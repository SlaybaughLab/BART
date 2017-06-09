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
(std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
  for (typename DoFHandler<dim>::active_cell_iterator
       cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
  {
    if (cell->is_locally_owned())
    {
      fv.reinit(cell);
      cell->get_dof_indices (local_dof_indices);
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            collision_at_qp[qi](i,j) = (fv->shape_value(i,qi) *
                                        fv->shape_value(j,qi));
      
      for (unsigned int qi=0; qi<n_q; ++qi)
        for (unsigned int i_dir=0; i_dir<n_dir; ++i_dir)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              streaming_at_qp[qi][i_dir](i,j) = ((fv->shape_grad(i,qi) *
                                                  omega_i[i_dir])
                                                 *
                                                 (fv->shape_grad(j,qi) *
                                                  omega_i[i_dir]));
      break;
    }// local cell
  }// cell
}

template <int dim>
void EvenParity<dim>::integrate_cell_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
  fv->reinit (cell);
  unsigned int material_id = cell->material_id ();
  for (unsigned int qi=0; qi<n_q; ++qi)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        cell_matrix(i,j) += (streaming_at_qp[qi][i_dir](i,j) *
                             all_inv_sigt[material_id][g]
                             +
                             collision_at_qp[qi](i,j) *
                             all_sigt[material_id][g]) * fv->JxW(qi);
      }
}

template <int dim>
void EvenParity<dim>::integrate_boundary_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 unsigned int &i_dir,
 unsigned int &g)
{
  for (unsigned int fn=0; fn<GeometryInfo<dim>::faces_per_cell; ++fn)
  {
    if (cell->at_boundary(fn) &&
        !is_reflective_bc[cell->face(fn)->boundary_id()])
    {
      fvf->reinit (cell,fn);
      const Tensor<1,dim> vec_n = fvf->normal_vector(0);
      double absndo = std::fabs (vec_n * omega_i[i_dir]);
      for (unsigned int qi=0; qi<n_qf; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (absndo *
                                 fvf->shape_value(i,qi) *
                                 fvf->shape_value(j,qi) *
                                 fvf->JxW(qi));
    }// non-ref bd
    else if (cell->at_boundary(fn) &&
             is_reflective_bc[cell->face(fn)->boundary_id()] &&
             is_explicit_reflective)
    {
      fvf.reinit (cell,fn);
      unsigned int boundary_id = cell->face(fn)->boundary_id ();
      unsigned int r_dir = get_reflective_direction_index (boundary_id, i_dir);
      const Tensor<1, dim> vec_n = fvf->normal_vector (0);
      double ndo = omega_i[i_dir] * vec_n;
      for (unsigned int qi=0; qi<n_qf; ++qi)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            cell_matrix(i,j) += (ndo *
                                 fvf->shape_value(i,qi) *
                                 (omega_i[r_dir] * fvf->shape_grad(j,qi)) *
                                 all_inv_sigt[material_id][g] *
                                 fvf->JxW(qi));
    }// explicit ref bd face
  }
}

template <int dim>
void EvenParity<dim>::integrate_interface_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,
 unsigned int &i_dir,
 unsigned int &g,
 FullMatrix<double> &vp_up,
 FullMatrix<double> &vp_un,
 FullMatrix<double> &vn_up,
 FullMatrix<double> &vn_un)
{
  fvf->reinit (cell,fn);
  fvf_nei->reinit (neig, cell->neighbor_face_no(fn));
  typename DoFHandler<dim>::cell_iterator neig = cell->neighbor(fn);
  
  const Tensor<1,dim> vec_n = fvf.normal_vector (0);
  
  double local_sigt = all_sigt[cell->material_id()][g];
  double local_inv_sigt = all_inv_sigt[cell->material_id()][g];
  double local_measure = cell->measure ();
  
  double neigh_sigt = all_sigt[neig->material_id()][g];
  double neigh_inv_sigt = all_inv_sigt[neig->material_id()][g];
  double neigh_measure = neig->measure ();
  
  double face_measure = cell->face(fn)->measure ();
  
  double avg_mfp_inv = 0.5 * (face_measure / (local_sigt * local_measure)
                              + face_measure / (neigh_sigt * neigh_measure));
  double sige = std::max(2.5, tensor_norms[i_dir] * c_penalty * avg_mfp_inv);
  
  for (unsigned int qi=0; qi<n_qf; ++qi)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
      {
        double theta = 1.0;
        vp_up(i,j) += (sige *
                       fvf->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       + theta *
                       (-
                        local_inv_sigt * vec_n * omega_i[i_dir] *
                        (omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                        fvf->shape_value(j,qi) * 0.5
                        -
                        local_inv_sigt * vec_n * omega_i[i_dir] *
                        fvf->shape_value(i,qi) * 0.5 *
                        (omega_i[i_dir] * fvf->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vp_un(i,j) += (-sige *
                       fvf->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       + theta *
                       (+
                        local_inv_sigt * vec_n * omega_i[i_dir] *
                        (omega_i[i_dir] * fvf->shape_grad(i,qi)) *
                        fvf_nei->shape_value(j,qi) * 0.5
                        -
                        neigh_inv_sigt * vec_n * omega_i[i_dir] *
                        fvf->shape_value(i,qi) * 0.5 *
                        (omega_i[i_dir] * fvf_nei->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vn_up(i,j) += (-sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf->shape_value(j,qi)
                       + theta *
                       (-
                        neigh_inv_sigt * vec_n * omega_i[i_dir] *
                        (omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                        fvf->shape_value(j,qi) * 0.5
                        +
                        local_inv_sigt * vec_n * omega_i[i_dir] *
                        fvf_nei->shape_value(i,qi) * 0.5 *
                        (omega_i[i_dir] * fvf->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
        
        vn_un(i,j) += (sige *
                       fvf_nei->shape_value(i,qi) *
                       fvf_nei->shape_value(j,qi)
                       + theta *
                       (+
                        neigh_inv_sigt * vec_n * omega_i[i_dir] *
                        (omega_i[i_dir] * fvf_nei->shape_grad(i,qi)) *
                        fvf_nei->shape_value(j,qi) * 0.5
                        +
                        neigh_inv_sigt * vec_n * omega_i[i_dir] *
                        fvf_nei->shape_value(i,qi) * 0.5 *
                        (omega_i[i_dir] * fvf_nei->shape_grad(j,qi)))
                       ) * fvf->JxW(qi);
      }

}

template class EvenParity<2>;
template class EvenParity<3>;
