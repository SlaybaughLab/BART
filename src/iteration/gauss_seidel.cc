#include "gauss_seidel.h"

template <int dim>
GaussSeidel<dim>::GaussSeidel (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    MGBase<dim> (prm, dat_ptr) {}

template <int dim>
GaussSeidel<dim>::~GaussSeidel () {}

template <int dim>
void GaussSeidel<dim>::NonthermalSolves (
    std::unique_ptr<EquationBase<dim>> &equ_ptr) {
  // loop over all nonthermal groups, assemble and solve one by one
  for (int g=0; g<this->g_thermal_; ++g) {
    equ_ptr->AssembleLinearForms (g);
    this->ig_ptr_->IGIterations (equ_ptr, g);
  }
}

template <int dim>
void GaussSeidel<dim>::ThermalIterations (
    std::unique_ptr<EquationBase<dim>> &equ_ptr) {
  const std::string equ_name = equ_ptr->GetEquName();
  double err = 1.0;
  while (err>this->err_phi_tol_) {
    double m = 0.0;
    for (int g=this->g_thermal_; g<this->n_group_; ++g) {
      // update scalar flux from previous mg iteration
      this->moments_prev_[std::make_tuple(g,0,0)] =
          this->mat_vec_->moments[equ_name][std::make_tuple(g,0,0)];
      // assemble rhs for current mg iteration at current group
      equ_ptr->AssembleLinearForms (g);
      // solve for current group
      this->ig_ptr_->IGIterations (equ_ptr, g);
      // calculate the error up to this group
      m = std::max (m,
          this->EstimatePhiDiff(this->mat_vec_->moments[equ_name][std::make_tuple(g,0,0)],
          this->moments_prev_[std::make_tuple(g,0,0)]));
    }
    err = m;
  }
}

template class GaussSeidel<1>;
template class GaussSeidel<2>;
template class GaussSeidel<3>;
