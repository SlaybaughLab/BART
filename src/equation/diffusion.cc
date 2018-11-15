#include "diffusion.h"

template <int dim>
Diffusion<dim>::Diffusion(const std::string &equation_name,
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    EquationBase<dim>(equation_name, prm, dat_ptr) {}

template <int dim>
Diffusion<dim>::~Diffusion () {}

template <int dim>
void Diffusion<dim>::PreassembleCellMatrices() {
  PreassembleCellMatricesDiffusion();
  PreassembleClosure();
}

template <int dim>
void Diffusion<dim>::PreassembleCellMatricesDiffusion () {
  auto cell = this->dat_ptr_->local_cells[0];
  this->fv_->reinit (cell);

  for (int qi=0; qi<this->n_q_; ++qi) {
    this->pre_collision_[qi] = dealii::FullMatrix<double> (
        this->dofs_per_cell_, this->dofs_per_cell_);
    pre_streaming_[qi] = dealii::FullMatrix<double> (
        this->dofs_per_cell_, this->dofs_per_cell_);
    for (int i=0; i<this->dofs_per_cell_; ++i)
      for (int j=0; j<this->dofs_per_cell_; ++j) {
        this->pre_collision_[qi](i,j) =
            this->fv_->shape_value(i,qi) *
            this->fv_->shape_value(j,qi);
        pre_streaming_[qi](i,j) =
            this->fv_->shape_grad(i,qi) *
            this->fv_->shape_grad(j,qi);
      }
  }
}

template <int dim>
void Diffusion<dim>::PreassembleClosure() {
  // Diffusion has zero closure for angular approximation
}

template <int dim>
void Diffusion<dim>::IntegrateCellBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::FullMatrix<double> &cell_matrix,
    const int &g,
    const int &dir) {
  int mid = cell->material_id ();
  auto siga = this->xsec_->sigt.at(mid)[g] - this->xsec_->sigs.at(mid)(g,g);
  for (int qi=0; qi<this->n_q_; ++qi)
    for (int i=0; i<this->dofs_per_cell_; ++i)
      for (int j=0; j<this->dofs_per_cell_; ++j)
        cell_matrix(i,j) += (pre_streaming_[qi](i,j) *
                             this->xsec_->diff_coef.at(mid)[g]
                             +
                             this->pre_collision_[qi](i,j) *
                             siga) * this->fv_->JxW(qi);
}

template <int dim>
void Diffusion<dim>::IntegrateScatteringLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q_);
  for (int gin=0; gin<this->n_group_; ++gin) {
    std::vector<double> local_flx (this->n_q_);
    this->fv_->get_function_values (
        this->mat_vec_->moments[this->equ_name_][std::make_tuple(gin, 0, 0)],
        local_flx);
    if (g!=gin && this->xsec_->sigs.at(mid)(gin,g)>bconst::kSmall)
      for (int qi=0; qi<this->n_q_; ++qi)
        q_at_qp[qi] += this->xsec_->sigs.at(mid)(gin,g) * local_flx[qi];
  }

  for (int qi=0; qi<this->n_q_; ++qi) {
    q_at_qp[qi] *= this->fv_->JxW(qi);
    for (int i=0; i<this->dofs_per_cell_; ++i)
      cell_rhs(i) += this->fv_->shape_value(i,qi) * q_at_qp[qi];
  }
}

template <int dim>
void Diffusion<dim>::IntegrateCellFixedLinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  int mid = cell->material_id ();
  std::vector<double> q_at_qp (this->n_q_);

  // retrieving "fixed" source at quadrature points
  if (!this->is_eigen_problem_) {
    q_at_qp = std::vector<double> (this->n_q_,
        this->xsec_->q.at(mid)[g]);
  } else if (this->is_eigen_problem_ && this->xsec_->is_material_fissile.at(mid)) {
    for (int gin=0; gin<this->n_group_; ++gin) {
      std::vector<double> local_flx (this->n_q_);
      this->fv_->get_function_values (
          this->mat_vec_->moments[this->equ_name_][std::make_tuple(gin, 0, 0)],
          local_flx);
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

template <int dim>
void Diffusion<dim>::IntegrateBoundaryBilinearForm (
    typename dealii::DoFHandler<dim>::active_cell_iterator &cell,
    const int &fn,/*face number*/
    dealii::Vector<double> &cell_rhs,
    const int &g,
    const int &dir) {
  int bd_id = cell->face(fn)->boundary_id ();
  const dealii::Tensor<1,dim> vec_n = this->fvf_->normal_vector(0);
  if (!this->have_reflective_bc_ || this->is_reflective_bc_.at(bd_id)) {
    // incident/vacuum boundary
    // TODO: check correctness
    for (int qi=0; qi<this->n_qf_; ++qi) {
      auto jxw = 0.25*this->fvf_->JxW(qi);
      for (int i=0; i<this->dofs_per_cell_; ++i)
        for (int j=0; j<this->dofs_per_cell_; ++j)
          cell_matrix(i,j) += (this->fvf_->shape_value(i,qi) *
                               this->fvf_->shape_value(j,qi) *
                               jxw);
    }
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;
