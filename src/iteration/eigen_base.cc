#include "eigen_base.h"

template <int dim>
EigenBase<dim>::EigenBase (const ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> dat_ptr)
    :
    IterationBase<dim>(prm, dat_ptr),
    err_k_tol_(1.0e-6),
    err_phi_tol_(1.0e-5) {}

template <int dim>
EigenBase<dim>::~EigenBase () {}

template <int dim>
void EigenBase<dim>::DoIterations (std::unordered_map<std::string,
    std::unique_ptr<EquationBase<dim>>> &equ_ptrs) {
  if (!this->do_nda_) {
    // assemble system matrices for transport equation
    equ_ptrs[this->ho_equ_name_]->AssembleBilinearForms ();
    // initialize fission process
    fiss_src_ = equ_ptrs[this->ho_equ_name_]->EstimateFissSrc ();
    keff_ = 1.0;
    // perform eigenvalue iterations
    EigenIterations (equ_ptrs[this->ho_equ_name_]);
  } else {
    // TODO: fill in NDA iterations
  }
}

// Override this function to do specific eigenvalue iteration as desired
template <int dim>
void EigenBase<dim>::EigenIterations (
    std::unique_ptr<EquationBase<dim>> &equ_ptr) {}

template <int dim>
void EigenBase<dim>::UpdatePrevSflxesFissSrcKeff (const std::string &equ_name) {
  // update scalar fluxes from previous eigen iteration
  for (unsigned int g=0; g<this->n_group_; ++g)
    this->moments_prev_[std::make_tuple(g,0,0)] = this->mat_vec_->moments[equ_name][std::make_tuple(g,0,0)];
  // update fission source from previous eigen iteration
  fiss_src_prev_ = fiss_src_;
  // update keff from previous eigen iteration
  keff_prev_ = keff_;
}

template <int dim>
void EigenBase<dim>::CalculateFissSrcKeff
(std::vector<Vector<double> > &sflxes_proc,
 std_cxx11::shared_ptr<EquationBase<dim> > equ_ptr)
{
  fiss_src = equ_ptr->estimate_fiss_src (sflxes_proc);
  keff = EstimateKeff ();
}

template <int dim>
inline double EigenBase<dim>::EstimateKeff () {
  return keff_prev_ * fiss_src_ / fiss_src_prev_;
}

template <int dim>
inline double EigenBase<dim>::EstimateKDiff () {
  return std::fabs (keff_ - keff_prev_) / keff_;
}

template <int dim>
double EigenBase<dim>::GetKeff () const {
  return keff_;
}

template class EigenBase<1>;
template class EigenBase<2>;
template class EigenBase<3>;
