#include "even_parity.h"

template <int dim>
EvenParity<dim>::EvenParity(const std::string &equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    EquationBase<dim>(equation_name, prm, dat_ptr) {
  if (this->discretization_=="dfem") {
    for (int i=0; i<this->n_dir_; ++i) {
      dealii::Tensor<2, dim> tmp =
          outer_product(this->omega_[i], this->omega_[i]);
      tensor_norms_.push_back(tmp.norm());
    }
    // note that the first constant on rhs is tunable and should be cautious
    // when tuning it.
    c_penalty_ = 1.0 * this->p_order_ * (this->p_order_ + 1.0);
  }
}

template <int dim>
EvenParity<dim>::~EvenParity () {}

template <int dim>
void EvenParity<dim>::PreassembleCellMatrices () {
  this->fv_->reinit (this->dat_ptr_->dof_handler.begin_active());

  for (int qi=0; qi<this->n_q_; ++qi) {
    this->pre_collision_[qi] =
        dealii::FullMatrix<double> (this->dofs_per_cell_, this->dofs_per_cell_);
    for (int i=0; i<this->dofs_per_cell_; ++i)
      for (int j=0; j<this->dofs_per_cell_; ++j)
        this->pre_collision_[qi](i,j) = (this->fv_->shape_value(i,qi) *
            this->fv_->shape_value(j,qi));
  }

  for (int qi=0; qi<this->n_q_; ++qi)
    for (int dir=0; dir<this->n_dir_; ++dir) {
      this->pre_streaming_[{dir, qi}] =
          dealii::FullMatrix<double> (this->dofs_per_cell_, this->dofs_per_cell_);
      for (int i=0; i<this->dofs_per_cell_; ++i)
        for (int j=0; j<this->dofs_per_cell_; ++j)
          this->pre_streaming_[{dir, qi}](i,j) =
             (this->fv_->shape_grad(i,qi) * this->omega_[dir])
             *
             (this->fv_->shape_grad(j,qi) * this->omega_[dir]);
    }
}

template <int dim>
void EvenParity<dim>::IntegrateCellBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::FullMatrix<double> &cell_matrix,
    const int &g,
    const int &dir) {
  int mid = cell->material_id ();
  for (int qi=0; qi<this->n_q_; ++qi)
    for (int i=0; i<this->dofs_per_cell_; ++i)
      for (int j=0; j<this->dofs_per_cell_; ++j)
        cell_matrix(i,j) += (this->pre_streaming_[{dir, qi}](i,j) *
                             this->xsec_->inv_sigt.at(mid)[g]
                             +
                             this->pre_collision_[qi](i,j) *
                             this->xsec_->sigt.at(mid)[g]) * this->fv_->JxW(qi);
}

template <int dim>
void EvenParity<dim>::IntegrateBoundaryBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    const int &fn,/*face number*/
    dealii::FullMatrix<double> &cell_matrix,
    const int &g,
    const int &dir) {
  int bd_id = cell->face(fn)->boundary_id ();
  const dealii::Tensor<1,dim> vec_n = this->fvf_->normal_vector(0);
  if (this->have_reflective_bc_ && this->is_reflective_bc_.at(bd_id)) {
    double inv_sigt =
        this->xsec_->inv_sigt.at(cell->material_id())[g];
    // hard coded part
    dealii::Tensor<1, dim> ref_angle =
        this->omega_[dir] - 2.0 * (this->omega_[dir] * vec_n) * vec_n;
    // standard part: not stable
    double ndo_inv_sigt = vec_n * this->omega_[dir] * inv_sigt;
    for (int qi=0; qi<this->n_qf_; ++qi)
      for (int i=0; i<this->dofs_per_cell_; ++i)
        for (int j=0; j<this->dofs_per_cell_; ++j)
          cell_matrix(i,j) += (- ndo_inv_sigt *
              this->fvf_->shape_value(i,qi) *
              (ref_angle * this->fvf_->shape_grad(j,qi)) *
              this->fvf_->JxW(qi));
  } else {
    // incident/vacuum boundary
    double absndo = std::fabs (vec_n * this->omega_[dir]);
    for (int qi=0; qi<this->n_qf_; ++qi)
      for (int i=0; i<this->dofs_per_cell_; ++i)
        for (int j=0; j<this->dofs_per_cell_; ++j)
          cell_matrix(i,j) += (absndo *
                               this->fvf_->shape_value(i,qi) *
                               this->fvf_->shape_value(j,qi) *
                               this->fvf_->JxW(qi));
  }
}

template <int dim>
void EvenParity<dim>::IntegrateBoundaryLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    const int &fn,/*face number*/
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  // We implement nothing for even parity boundary linear form. In reflective BC,
  // even parity realizes it via modifying bilinear form. Also, we assume vacuum
  // BC so nothing needs to be realized therein.
}

template <int dim>
void EvenParity<dim>::IntegrateInterfaceBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    typename dealii::DoFHandler<dim>::cell_iterator &neigh,/*cell iterator for cell*/
    const int &fn,/*concerning face number in local cell*/
    dealii::FullMatrix<double> &vp_up,
    dealii::FullMatrix<double> &vp_un,
    dealii::FullMatrix<double> &vn_up,
    dealii::FullMatrix<double> &vn_un,
    const int &g,
    const int &dir) {
  // TODO: vec_n below only works if the mesh is linear as normal vectors point
  // to the same direction per cell face in such a case.
  const dealii::Tensor<1,dim> vec_n = this->fvf_->normal_vector (0);
  int mid = cell->material_id ();
  int mid_nei = neigh->material_id ();
  double local_sigt = this->xsec_->sigt.at(mid)[g];
  double local_inv_sigt = this->xsec_->inv_sigt.at(mid)[g];
  double neigh_sigt = this->xsec_->sigt.at(mid_nei)[g];
  double neigh_inv_sigt = this->xsec_->inv_sigt.at(mid_nei)[g];
  /* The following section would not work with 1D
  double local_measure = cell->measure ();
  double neigh_measure = neigh->measure ();
  double face_measure = cell->face(fn)->measure ();

  double avg_mfp_inv = 0.5 * (face_measure / (local_sigt * local_measure)
      + face_measure / (neigh_sigt * neigh_measure));
  */
  double local_measure = cell->diameter ();
  double neigh_measure = neigh->diameter ();
  double avg_mfp_inv = 0.5 / (local_sigt * local_measure)
      + 0.5 / (neigh_sigt * neigh_measure);
  // the min value of sige is another tricky part and one should use cautions
  double sige = std::max(0.25, tensor_norms_[dir] * c_penalty_ * avg_mfp_inv);

  double half_ndo = 0.5 * vec_n * this->omega_[dir];
  //double sige = std::max(std::fabs (ndo),0.25);
  for (int qi=0; qi<this->n_qf_; ++qi)
    for (int i=0; i<this->dofs_per_cell_; ++i)
      for (int j=0; j<this->dofs_per_cell_; ++j) {
        vp_up(i,j) += (sige *
            this->fvf_->shape_value(i,qi) *
            this->fvf_->shape_value(j,qi)
            -
            local_inv_sigt * half_ndo *
            (this->omega_[dir] * this->fvf_->shape_grad(i,qi)) *
            this->fvf_->shape_value(j,qi)
            -
            local_inv_sigt * half_ndo *
            this->fvf_->shape_value(i,qi) *
            (this->omega_[dir] * this->fvf_->shape_grad(j,qi)))
            * this->fvf_->JxW(qi);

        vp_un(i,j) += (-sige *
            this->fvf_->shape_value(i,qi) *
            this->fvf_nei_->shape_value(j,qi)
            +
            local_inv_sigt * half_ndo *
            (this->omega_[dir] * this->fvf_->shape_grad(i,qi)) *
            this->fvf_nei_->shape_value(j,qi)
            -
            neigh_inv_sigt * half_ndo *
            this->fvf_->shape_value(i,qi) *
            (this->omega_[dir] * this->fvf_nei_->shape_grad(j,qi))
            ) * this->fvf_->JxW(qi);

        vn_up(i,j) += (-sige *
            this->fvf_nei_->shape_value(i,qi) *
            this->fvf_->shape_value(j,qi)
            -
            neigh_inv_sigt * half_ndo *
            (this->omega_[dir] * this->fvf_nei_->shape_grad(i,qi)) *
            this->fvf_->shape_value(j,qi)
            +
            local_inv_sigt * half_ndo *
            this->fvf_nei_->shape_value(i,qi) *
            (this->omega_[dir] * this->fvf_->shape_grad(j,qi))
            ) * this->fvf_->JxW(qi);

        vn_un(i,j) += (sige *
            this->fvf_nei_->shape_value(i,qi) *
            this->fvf_nei_->shape_value(j,qi)
            +
            neigh_inv_sigt * half_ndo *
            (this->omega_[dir] * this->fvf_nei_->shape_grad(i,qi)) *
            this->fvf_nei_->shape_value(j,qi)
            +
            neigh_inv_sigt * half_ndo *
            this->fvf_nei_->shape_value(i,qi) *
            (this->omega_[dir] * this->fvf_nei_->shape_grad(j,qi))
            ) * this->fvf_->JxW(qi);
      }

}

template <int dim>
void EvenParity<dim>::IntegrateScatteringLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  // dir info is irrelavant for even parity
  // TODO: Anisotropic scattering is not included
  int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q_);
  for (int gin=0; gin<this->n_group_; ++gin) {
    std::vector<double> local_flx (this->n_q_);
    this->fv_->get_function_values (
        this->mat_vec_->moments[this->equ_name_][{gin, 0, 0}], local_flx);
    for (int qi=0; qi<this->n_q_; ++qi)
      q_at_qp[qi] += (this->xsec_->sigs_per_ster.at(mid)(gin,g) *
                      local_flx[qi]);
  }

  for (int qi=0; qi<this->n_q_; ++qi) {
    q_at_qp[qi] *= this->fv_->JxW(qi);
    for (int i=0; i<this->dofs_per_cell_; ++i)
      cell_rhs(i) += this->fv_->shape_value(i,qi) * q_at_qp[qi];
  }
}

template <int dim>
void EvenParity<dim>::IntegrateCellFixedLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  // dir info is irrelavant for even parity
  int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q_);

  // retrieving "fixed" source at quadrature points
  if (!this->is_eigen_problem_) {
    q_at_qp = std::vector<double> (this->n_q_,
        this->xsec_->q_per_ster.at(mid)[g]);
  } else if (this->is_eigen_problem_ && this->xsec_->is_material_fissile.at(mid)) {
    for (int gin=0; gin<this->n_group_; ++gin)
    {
      std::vector<double> local_flx (this->n_q_);
      this->fv_->get_function_values (
          this->mat_vec_->moments[this->equ_name_][{gin, 0, 0}], local_flx);
      for (int qi=0; qi<this->n_q_; ++qi)
        q_at_qp[qi] += this->scaled_fiss_transfer_.at(mid)(gin, g) *
            local_flx[qi];
    }
  }
  // calculate cell rhs:
  for (int qi=0; qi<this->n_q_; ++qi) {
    q_at_qp[qi] *= this->fv_->JxW(qi);
    for (int i=0; i<this->dofs_per_cell_; ++i)
      cell_rhs(i) += this->fv_->shape_value(i,qi) * q_at_qp[qi];
  }
}

template class EvenParity<1>;
template class EvenParity<2>;
template class EvenParity<3>;
