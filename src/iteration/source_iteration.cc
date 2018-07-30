#include "source_iteration.h"

template <int dim>
SourceIteration<dim>::SourceIteration (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    IGBase<dim> (prm, dat_ptr) {}

template <int dim>
SourceIteration<dim>::~SourceIteration () {}

template <int dim>
void SourceIteration<dim>::IGIterations (
    std::unique_ptr<EquationBase<dim>> &equ_ptr,
    const int &g) {
  double err = 1.0;
  int iter = 0;
  const std::string equ_name = equ_ptr->GetEquName();
  while (err>this->err_phi_tol_) {
    // generate rhs for group g
    equ_ptr->AssembleLinearForms (g);
    // solve all the directions in group g
    equ_ptr->SolveInGroup (g);
    // generate moments
    equ_ptr->GenerateMoments (this->mat_vec_->moments[equ_name],
        this->moments_prev_, g);
    // calculate the difference of moments for convergence check
    err = this->EstimatePhiDiff (this->mat_vec_->moments[equ_name],
        this->moments_prev_);
    //this->pcout << "SI error: " << err << ", Iter: " << iter++ << std::endl;
  }
}

template class SourceIteration<1>;
template class SourceIteration<2>;
template class SourceIteration<3>;
