#include "power_iteration.h"

template <int dim>
PowerIteration<dim>::PowerIteration (const ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> dat_ptr)
    :
    EigenBase<dim> (prm, dat_ptr) {}

template <int dim>
PowerIteration<dim>::~PowerIteration () {}

template <int dim>
void PowerIteration<dim>::EigenIterations (
    std::unique_ptr<EquationBase<dim>> equ_ptr) {
  // TODO: nda requires overriding do_iterations
  const std::string equ_name = equ_ptr->GetEquName();
  double err_k = 1.0;
  double err_phi = 1.0;
  int iter = 0;
  while (err_k>this->err_k_tol_ || err_phi>this->err_phi_tol_) {
    // update sflxes, fission source and keff from previous fission with current sflxes
    this->UpdatePrevSflxesFissSrcKeff (equ_name);
    // scale chi_nu_sigf by a factor of keff
    equ_ptr->ScaleFissTransferMatrices (this->keff_);
    // assemble fission source as a "fixed source"
    equ_ptr->AssembleFixedLinearForms ();
    // perform multigroup iterations
    this->mg_ptr_->MGIterations (equ_ptr);
    // calculate fission source and keff thereafter
    this->CalculateFissSrcKeff (equ_ptr);
    // calculate errors of QoI for convergence check
    err_phi = this->EstimatePhiDiff (this->mat_vec_->moments[equ_name],
        this->moments_prev_);
    err_k = this->EstimateKDiff ();
    // print on screen about the errors
    //this->pcout << std::endl << std::endl << "PI iter: " << iter++
    //<< ", err_k: " << err_k << ", err_phi: " << err_phi << std::endl << std::endl;
  }
}

template class PowerIteration<1>;
template class PowerIteration<2>;
template class PowerIteration<3>;
