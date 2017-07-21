#include "transport_base.h"
#include "even_parity.h"

template <int dim>
EvenParity<dim>::EvenParity (ParameterHandler &prm)
:
TransportBase<dim>(prm)
{
}

template <int dim>
EvenParity<dim>::~EvenParity ()
{
}

template <int dim>
void EvenParity<dim>::pre_assemble_cell_matrices
(const std_cxx11::shared_ptr<FEValues<dim> > fv,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
        collision_at_qp[qi](i,j) = (fv->shape_value(i,qi) *
                                    fv->shape_value(j,qi));

  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i_dir=0; i_dir<this->n_dir; ++i_dir)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          streaming_at_qp[qi][i_dir](i,j) = ((fv->shape_grad(i,qi) *
                                              this->omega_i[i_dir])
                                             *
                                             (fv->shape_grad(j,qi) *
                                              this->omega_i[i_dir]));
}

template <int dim>
void EvenParity<dim>::integrate_cell_bilinear_form
(const std_cxx11::shared_ptr<FEValues<dim> > fv,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
  unsigned int mid = cell->material_id ();
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
        cell_matrix(i,j) += (streaming_at_qp[qi][i_dir](i,j) *
                             this->all_inv_sigt[mid][g]
                             +
                             collision_at_qp[qi](i,j) *
                             this->all_sigt[mid][g]) * fv->JxW(qi);
}

template <int dim>
void EvenParity<dim>::integrate_boundary_bilinear_form
(const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g)
{
  unsigned int bd_id = cell->face(fn)->boundary_id ();
  const Tensor<1,dim> vec_n = fvf->normal_vector(0);
  if (this->have_reflective_bc && this->is_reflective_bc[bd_id])
  {
    unsigned int inv_sigt = this->all_inv_sigt[cell->material_id()][g];
    // hard coded part
    Tensor<1, dim> ref_angle =
    this->omega_i[i_dir] - 2.0 * (this->omega_i[i_dir] * vec_n) * vec_n;
    // standard part: not stable
    //unsigned int r_dir = this->get_reflective_direction_index (bd_id, i_dir);
    //ref_angle = this->omega_i[r_dir];
    double ndo_inv_sigt = vec_n * this->omega_i[i_dir] * inv_sigt;
    for (unsigned int qi=0; qi<this->n_qf; ++qi)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          cell_matrix(i,j) += (- ndo_inv_sigt *
                               fvf->shape_value(i,qi) *
                               (ref_angle * fvf->shape_grad(j,qi)) *
                               fvf->JxW(qi));
  }
  else/* if (!this->have_reflective_bc ||
           (this->have_reflective_bc && !this->is_reflective_bc[bd_id]))*/
  {
    double absndo = std::fabs (vec_n * this->omega_i[i_dir]);
    for (unsigned int qi=0; qi<this->n_qf; ++qi)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          cell_matrix(i,j) += (absndo *
                               fvf->shape_value(i,qi) *
                               fvf->shape_value(j,qi) *
                               fvf->JxW(qi));
  }// non-ref bd
}

template <int dim>
void EvenParity<dim>::integrate_interface_bilinear_form
(const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf,
 const std_cxx11::shared_ptr<FEFaceValues<dim> > fvf_nei,
 typename DoFHandler<dim>::active_cell_iterator &cell,
 typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
 unsigned int &fn,/*concerning face number in local cell*/
 unsigned int &i_dir,
 unsigned int &g,
 FullMatrix<double> &vp_up,
 FullMatrix<double> &vp_un,
 FullMatrix<double> &vn_up,
 FullMatrix<double> &vn_un)
{
  const Tensor<1,dim> vec_n = fvf->normal_vector (0);
  unsigned int mid = cell->material_id ();
  unsigned int mid_nei = neigh->material_id ();
  double local_sigt = this->all_sigt[mid][g];
  double local_inv_sigt = this->all_inv_sigt[mid][g];
  double local_measure = cell->measure ();
  double neigh_sigt = this->all_sigt[mid_nei][g];
  double neigh_inv_sigt = this->all_inv_sigt[mid_nei][g];
  double neigh_measure = neigh->measure ();
  double face_measure = cell->face(fn)->measure ();

  double avg_mfp_inv = 0.5 * (face_measure / (local_sigt * local_measure)
                              + face_measure / (neigh_sigt * neigh_measure));
  double sige = std::max(0.25, this->tensor_norms[i_dir] * this->c_penalty * avg_mfp_inv);

  double half_ndo = 0.5 * vec_n * this->omega_i[i_dir];
  //double sige = std::max(std::fabs (ndo),0.25);
  for (unsigned int qi=0; qi<this->n_qf; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
      {
        vp_up(i,j) += (sige *
                       fvf->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       -
                       local_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                       fvf->shape_value(j,qi)
                       -
                       local_inv_sigt * half_ndo *
                       fvf->shape_value(i,qi) *
                       (this->omega_i[i_dir] * fvf->shape_grad(j,qi))
                       ) * fvf->JxW(qi);

        vp_un(i,j) += (-sige *
                       fvf->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       +
                       local_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                       fvf_nei->shape_value(j,qi)
                       -
                       neigh_inv_sigt * half_ndo *
                       fvf->shape_value(i,qi) *
                       (this->omega_i[i_dir] * fvf_nei->shape_grad(j,qi))
                       ) * fvf->JxW(qi);

        vn_up(i,j) += (-sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       -
                       neigh_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                       fvf->shape_value(j,qi)
                       +
                       local_inv_sigt * half_ndo *
                       fvf_nei->shape_value(i,qi) *
                       (this->omega_i[i_dir] * fvf->shape_grad(j,qi))
                       ) * fvf->JxW(qi);

        vn_un(i,j) += (sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       +
                       neigh_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                       fvf_nei->shape_value(j,qi)
                       +
                       neigh_inv_sigt * half_ndo *
                       fvf_nei->shape_value(i,qi) *
                       (this->omega_i[i_dir] * fvf_nei->shape_grad(j,qi))
                       ) * fvf->JxW(qi);
      }

}

template <int dim>
void EvenParity<dim>::generate_ho_rhs ()
{
  for (unsigned int g=0; g<this->n_group; ++g)
    for (unsigned int i_dir=0; i_dir<this->n_dir; ++i_dir)
    {
      unsigned int k = this->get_component_index (i_dir, g);
      if (i_dir==0 && !this->do_nda)
      {
        *(this->vec_ho_rhs[k]) = 0.0;
        for (unsigned int ic=0; ic<this->local_cells.size (); ++ic)
        {
          Vector<double> cell_rhs (this->dofs_per_cell);
          typename DoFHandler<dim>::active_cell_iterator cell = this->local_cells[ic];
          cell->get_dof_indices (this->local_dof_indices);
          this->fv->reinit (cell);
          unsigned int mid = cell->material_id ();
          std::vector<std::vector<double> > local_sflxes
          (this->n_group, std::vector<double>(this->n_q));
          for (unsigned int gin=0; gin<this->n_group; ++gin)
            this->fv->get_function_values (this->sflx_proc[gin], local_sflxes[gin]);
          
          for (unsigned int qi=0; qi<this->n_q; ++qi)
          {
            double q_at_qp = 0.0;
            for (unsigned int gin=0; gin<this->n_group; ++gin)
              q_at_qp += (this->all_sigs_per_ster[mid][gin][g]<1.0e-13?0.0:
                          (this->all_sigs_per_ster[mid][gin][g] * local_sflxes[gin][qi]));
            for (unsigned int i=0; i<this->dofs_per_cell; ++i)
              cell_rhs (i) += this->vec_test_at_qp[ic](qi, i) * q_at_qp;
          }
          this->vec_ho_rhs[k]->add (this->local_dof_indices, cell_rhs);
        }// local cells
        this->vec_ho_rhs[k]->compress (VectorOperation::add);
        *(this->vec_ho_rhs[k]) += *(this->vec_ho_fixed_rhs[k]);
      }// zeroth direction per group
      else
        *(this->vec_ho_rhs[k]) = *(this->vec_ho_rhs[this->get_component_index(0, g)]);
    // Note that reflective boundary condition is carreid out using explicit reflective
    // algorithm. See Memo 2 for details.
    }// i_dir
}

template <int dim>
void EvenParity<dim>::generate_ho_fixed_source ()
{
  for (unsigned int g=0; g<this->n_group; ++g)
    for (unsigned int i_dir=0; i_dir<this->n_dir; ++i_dir)
    {
      unsigned int k = this->get_component_index (i_dir, g);
      if (i_dir==0)
      {
        *(this->vec_ho_fixed_rhs[k]) = 0.0;
        for (unsigned int ic=0; ic<this->local_cells.size (); ++ic)
        {
          Vector<double> cell_rhs (this->dofs_per_cell);
          typename DoFHandler<dim>::active_cell_iterator cell = this->local_cells[ic];
          unsigned int mid = cell->material_id ();
          
          if ((this->is_eigen_problem && this->is_material_fissile[mid]) ||
              (!this->is_eigen_problem &&
               (this->do_nda || (!this->do_nda && this->all_q_per_ster[mid][g]>1.0e-13))))
          {
            this->fv->reinit (cell);
            cell->get_dof_indices (this->local_dof_indices);
            std::vector<std::vector<double> > local_sflxes (this->n_group, std::vector<double>(this->n_q));
            for (unsigned int gin=0; gin<this->n_group; ++gin)
            {
              if (this->do_nda)
                this->fv->get_function_values (this->lo_sflx_proc[gin], local_sflxes[gin]);
              else if (!this->do_nda && this->is_eigen_problem)
                this->fv->get_function_values (this->sflx_proc_prev_gen[gin], local_sflxes[gin]);
            }
            
            for (unsigned int qi=0; qi<this->n_q; ++qi)
            {
              double q_at_qp = 0.0;
              // calculate pointwise source per spatial quadrature point
              if (this->do_nda)
              {
                if (this->is_eigen_problem)
                  for (unsigned int gin=0; gin<this->n_group; ++gin)
                    q_at_qp += (this->scat_scaled_fiss_transfer_per_ster[mid][gin][g]<1.0e-13?0.0:
                                (this->scat_scaled_fiss_transfer_per_ster[mid][gin][g] *
                                 local_sflxes[gin][qi]));
                else
                  for (unsigned int gin=0; gin<this->n_group; ++gin)
                    q_at_qp += (this->all_sigs_per_ster[mid][gin][g]<1.0e-13?0.0:
                                (this->all_sigs_per_ster[mid][gin][g] *
                                 local_sflxes[gin][qi]));
              }
              else// no NDA
              {
                if (this->is_eigen_problem)// fission source is the fixed source
                  for (unsigned int gin=0; gin<this->n_group; ++gin)
                    q_at_qp += (!this->is_material_fissile[mid]?0.0:
                                (this->scaled_fiss_transfer_per_ster[mid][gin][g] *
                                 local_sflxes[gin][qi]));
                else
                  q_at_qp += this->all_q_per_ster[mid][g];
              }
              for (unsigned int i=0; i<this->dofs_per_cell; ++i)
                cell_rhs (i) += this->vec_test_at_qp[ic](qi, i) * q_at_qp;
            }
            this->vec_ho_fixed_rhs[k]->add (this->local_dof_indices, cell_rhs);
          }// when to calculate rhs
        }// loop over local cells
        this->vec_ho_fixed_rhs[k]->compress (VectorOperation::add);
      }// first direction per group
      else
        *(this->vec_ho_fixed_rhs[k]) =
        *(this->vec_ho_fixed_rhs[this->get_component_index(0, g)]);
    }
}

template class EvenParity<2>;
template class EvenParity<3>;
