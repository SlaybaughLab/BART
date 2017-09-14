#include "even_parity.h"

template <int dim>
EvenParity<dim>::EvenParity
(std::string equation_name,
 const ParameterHandler &prm,
 const std_cxx11::shared_ptr<MeshGenerator<dim> > msh_ptr,
 const std_cxx11::shared_ptr<AQBase<dim> > aqd_ptr,
 const std_cxx11::shared_ptr<MaterialProperties> mat_ptr)
:
EquationBase<dim>(equation_name, prm, msh_ptr, aqd_ptr, mat_ptr)
{
  AssertThrow(equation_name=="ep",
              ExcMessage("equation built incorrectly"));
  if (this->discretization=="dfem")
  {
    for (unsigned int i=0; i<this->n_dir; ++i)
    {
      Tensor<2, dim> tensor_tmp = outer_product(this->omega_i[i], this->omega_i[i]);
      tensor_norms.push_back(tensor_tmp.norm());
    }
    // note that the first constant on rhs is tunable and should be cautious
    // when tuning it.
    c_penalty = 1.0 * this->p_order * (this->p_order + 1.0);
  }
}

template <int dim>
EvenParity<dim>::~EvenParity ()
{
}

template <int dim>
void EvenParity<dim>::pre_assemble_cell_matrices
(typename DoFHandler<dim>::active_cell_iterator &cell,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp)
{
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
        collision_at_qp[qi](i,j) = (this->fv->shape_value(i,qi) *
                                    this->fv->shape_value(j,qi));
  
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i_dir=0; i_dir<this->n_dir; ++i_dir)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          streaming_at_qp[qi][i_dir](i,j) = ((this->fv->shape_grad(i,qi) *
                                              this->omega_i[i_dir])
                                             *
                                             (this->fv->shape_grad(j,qi) *
                                              this->omega_i[i_dir]));
}

template <int dim>
void EvenParity<dim>::integrate_cell_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 FullMatrix<double> &cell_matrix,
 std::vector<std::vector<FullMatrix<double> > > &streaming_at_qp,
 std::vector<FullMatrix<double> > &collision_at_qp,
 const unsigned int &g,
 const unsigned int &i_dir)
{
  unsigned int mid = cell->material_id ();
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
        cell_matrix(i,j) += (streaming_at_qp[qi][i_dir](i,j) *
                             this->all_inv_sigt[mid][g]
                             +
                             collision_at_qp[qi](i,j) *
                             this->all_sigt[mid][g]) * this->fv->JxW(qi);
}

template <int dim>
void EvenParity<dim>::integrate_boundary_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 unsigned int &fn,/*face number*/
 FullMatrix<double> &cell_matrix,
 const unsigned int &g,
 const unsigned int &i_dir)
{
  unsigned int bd_id = cell->face(fn)->boundary_id ();
  const Tensor<1,dim> vec_n = this->fvf->normal_vector(0);
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
                               this->fvf->shape_value(i,qi) *
                               (ref_angle * this->fvf->shape_grad(j,qi)) *
                               this->fvf->JxW(qi));
  }
  else/* if (!this->have_reflective_bc ||
       (this->have_reflective_bc && !this->is_reflective_bc[bd_id]))*/
  {
    double absndo = std::fabs (vec_n * this->omega_i[i_dir]);
    for (unsigned int qi=0; qi<this->n_qf; ++qi)
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
        for (unsigned int j=0; j<this->dofs_per_cell; ++j)
          cell_matrix(i,j) += (absndo *
                               this->fvf->shape_value(i,qi) *
                               this->fvf->shape_value(j,qi) *
                               this->fvf->JxW(qi));
  }// non-ref bd
}

template <int dim>
void EvenParity<dim>::integrate_interface_bilinear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 typename DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
 unsigned int &fn,/*concerning face number in local cell*/
 FullMatrix<double> &vp_up,
 FullMatrix<double> &vp_un,
 FullMatrix<double> &vn_up,
 FullMatrix<double> &vn_un,
 const unsigned int &g,
 const unsigned int &i_dir)
{
  const Tensor<1,dim> vec_n = this->fvf->normal_vector (0);
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
  // the min value of sige is another tricky part and one should use cautions
  double sige = std::max(0.25, tensor_norms[i_dir] * c_penalty * avg_mfp_inv);
  
  double half_ndo = 0.5 * vec_n * this->omega_i[i_dir];
  //double sige = std::max(std::fabs (ndo),0.25);
  for (unsigned int qi=0; qi<this->n_qf; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
      {
        vp_up(i,j) += (sige *
                       this->fvf->shape_value(i,qi) *
                       this->fvf->shape_value(j,qi)
                       -
                       local_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * this->fvf->shape_grad(i,qi)) *
                       this->fvf->shape_value(j,qi)
                       -
                       local_inv_sigt * half_ndo *
                       this->fvf->shape_value(i,qi) *
                       (this->omega_i[i_dir] * this->fvf->shape_grad(j,qi))
                       ) * this->fvf->JxW(qi);
        
        vp_un(i,j) += (-sige *
                       this->fvf->shape_value(i,qi) *
                       this->fvf_nei->shape_value(j,qi)
                       +
                       local_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * this->fvf->shape_grad(i,qi)) *
                       this->fvf_nei->shape_value(j,qi)
                       -
                       neigh_inv_sigt * half_ndo *
                       this->fvf->shape_value(i,qi) *
                       (this->omega_i[i_dir] * this->fvf_nei->shape_grad(j,qi))
                       ) * this->fvf->JxW(qi);
        
        vn_up(i,j) += (-sige *
                       this->fvf_nei->shape_value(i,qi) *
                       this->fvf->shape_value(j,qi)
                       -
                       neigh_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * this->fvf_nei->shape_grad(i,qi)) *
                       this->fvf->shape_value(j,qi)
                       +
                       local_inv_sigt * half_ndo *
                       this->fvf_nei->shape_value(i,qi) *
                       (this->omega_i[i_dir] * this->fvf->shape_grad(j,qi))
                       ) * this->fvf->JxW(qi);
        
        vn_un(i,j) += (sige *
                       this->fvf_nei->shape_value(i,qi) *
                       this->fvf_nei->shape_value(j,qi)
                       +
                       neigh_inv_sigt * half_ndo *
                       (this->omega_i[i_dir] * this->fvf_nei->shape_grad(i,qi)) *
                       this->fvf_nei->shape_value(j,qi)
                       +
                       neigh_inv_sigt * half_ndo *
                       this->fvf_nei->shape_value(i,qi) *
                       (this->omega_i[i_dir] * this->fvf_nei->shape_grad(j,qi))
                       ) * this->fvf->JxW(qi);
      }
  
}

template <int dim>
void EvenParity<dim>::integrate_scattering_linear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 Vector<double> &cell_rhs,
 std::vector<Vector<double> > &sflx_proc,
 const unsigned int &g,
 const unsigned int &i_dir)
{
  // i_dir info is irrelavant for even parity
  unsigned int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q);
  for (unsigned int gin=0; gin<this->n_group; ++gin)
  {
    std::vector<double> local_flx (this->n_q);
    this->fv->get_function_values (sflx_proc[gin], local_flx);
    for (unsigned int qi=0; qi<this->n_q; ++qi)
      q_at_qp[qi] += (this->all_sigs_per_ster[mid][gin][g] *
                      local_flx[qi]);
  }

  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      cell_rhs(i) += (this->fv->shape_value(i,qi) *
                      q_at_qp[qi] *
                      this->fv->JxW(qi));
}

template <int dim>
void EvenParity<dim>::integrate_cell_fixed_linear_form
(typename DoFHandler<dim>::active_cell_iterator &cell,
 Vector<double> &cell_rhs,
 std::vector<Vector<double> > &sflx_prev,
 const unsigned int &g,
 const unsigned int &i_dir)
{
  // i_dir info is irrelavant for even parity
  unsigned int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q);
  
  // retrieving "fixed" source at quadrature points
  if (!this->is_eigen_problem)
    q_at_qp = std::vector<double> (this->n_q, this->all_q_per_ster[mid][g]);
  else if (this->is_eigen_problem && this->is_material_fissile[mid])
    for (unsigned int gin=0; gin<this->n_q; ++gin)
    {
      std::vector<double> local_flx (this->n_q);
      this->fv->get_function_values (sflx_prev[gin], local_flx);
      for (unsigned int qi=0; qi<this->n_q; ++qi)
        q_at_qp[qi] += (this->scaled_fiss_transfer_per_ster[mid][gin][g] *
                        local_flx[qi]);
    }
  
  // calculate cell rhs:
  for (unsigned int qi=0; qi<this->n_q; ++qi)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      cell_rhs(i) += (this->fv->shape_value(i,qi) *
                      q_at_qp[qi] *
                      this->fv->JxW(qi));
}

template class EvenParity<2>;
template class EvenParity<3>;
