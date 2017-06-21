#include "../include/transport_base.h"
#include "../include/even_parity.h"

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
  unsigned int boundary_id = cell->face(fn)->boundary_id ();
  
  if (!this->is_reflective_bc[boundary_id])
  {
    const Tensor<1,dim> vec_n = fvf->normal_vector(0);
    double absndo = std::fabs (vec_n * this->omega_i[i_dir]);
    for (unsigned int qi=0; qi<this->n_qf; ++qi)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          cell_matrix(i,j) += (absndo *
                               fvf->shape_value(i,qi) *
                               fvf->shape_value(j,qi) *
                               fvf->JxW(qi));
  }// non-ref bd
  else if (this->is_reflective_bc[boundary_id] &&
           this->is_explicit_reflective)
  {
    unsigned int mid = cell->material_id ();
    unsigned int r_dir = this->get_reflective_direction_index (boundary_id, i_dir);
    const Tensor<1, dim> vec_n = fvf->normal_vector (0);
    double ndo = this->omega_i[i_dir] * vec_n;
    for (unsigned int qi=0; qi<this->n_qf; ++qi)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          cell_matrix(i,j) += (ndo *
                               fvf->shape_value(i,qi) *
                               (this->omega_i[r_dir] * fvf->shape_grad(j,qi)) *
                               this->all_inv_sigt[mid][g] *
                               fvf->JxW(qi));
  }// explicit ref bd face
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
  
  double ndo = vec_n * this->omega_i[i_dir];
  //double sige = std::max(std::fabs (ndo),0.25);
  //std::cout << "sige: " << sige << std::endl;
  
  for (unsigned int qi=0; qi<this->n_qf; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
      {
        double theta = 1.0;
        vp_up(i,j) += (sige *
                       fvf->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       + theta*
                       (-
                        local_inv_sigt * ndo *
                        (this->omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                        fvf->shape_value(j,qi) * 0.5
                        -
                        local_inv_sigt * ndo *
                        fvf->shape_value(i,qi) * 0.5 *
                        (this->omega_i[i_dir] * fvf->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vp_un(i,j) += (-sige *
                       fvf->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       + theta*
                       (+
                        local_inv_sigt * ndo *
                        (this->omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                        fvf_nei->shape_value(j,qi) * 0.5
                        -
                        neigh_inv_sigt * ndo *
                        fvf->shape_value(i,qi) * 0.5 *
                        (this->omega_i[i_dir] * fvf_nei->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vn_up(i,j) += (-sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       + theta*
                       (-
                        neigh_inv_sigt * ndo *
                        (this->omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                        fvf->shape_value(j,qi) * 0.5
                        +
                        local_inv_sigt * ndo *
                        fvf_nei->shape_value(i,qi) * 0.5 *
                        (this->omega_i[i_dir] * fvf->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vn_un(i,j) += (sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       + theta*
                       (+
                        neigh_inv_sigt * ndo *
                        (this->omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                        fvf_nei->shape_value(j,qi) * 0.5
                        +
                        neigh_inv_sigt * ndo *
                        fvf_nei->shape_value(i,qi) * 0.5 *
                        (this->omega_i[i_dir] * fvf_nei->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
      }
  
}

template class EvenParity<2>;
template class EvenParity<3>;
