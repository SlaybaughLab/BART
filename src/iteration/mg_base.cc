#include "mg_base.h"

template <int dim>
MGBase<dim>::MGBase (const ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> dat_ptr)
    :
    IterationBase<dim> (prm, dat_ptr),
    g_thermal_(prm.get_integer("thermal group boundary")),
    err_phi_tol_(1.0e-5) {
  AssertThrow (g_thermal<this->n_group && g_thermal>=0,
               ExcMessage("Invalid thermal upper boundary"));

  // TODO: needs change if anisotropic scattering is needed
  sflxes_proc_prev_mg.resize (this->n_group);
}

template <int dim>
MGBase<dim>::~MGBase () {}

template <int dim>
void MGBase<dim>::DoIterations (std::unordered_map<std::string,
    std::unique_ptr<EquationBase<dim>>> &equ_ptrs) {
  // this function will be called only when it's not an eigenvalue problem. Copy
  // the following assertion in ALL derived class of MGBase.
  AssertThrow (!this->is_eigen_problem_,
               ExcMessage("This function shall not be called if it's eigen prob."));
  // purpose of overriding do_iteration is such that we provide another outer
  // iteration here if NDA is used
  if (this->do_nda_) {
    // TODO: fill in NDA iterations using MGIterations
  } else {
    // assemble bilinear forms of available equations
    equ_ptrs[this->ho_equ_name_]->AssembleBilinearForms ();
    // multigroup iterations. Note we would not need to override mg_iterations
    // until we want to do JFNK
    this->MGIterations (equ_ptrs[this->ho_equ_name_]);
  }
}

template <int dim>
void MGBase<dim>::MGIterations (std::unique_ptr<EquationBase<dim>> equ_ptr) {
  // this function needs to be overridden if JFNK is desired
  // solve for epithermal and fast groups
  NonthermalSolves (equ_ptr);
  // thermal group iterations
  ThermalIterations (equ_ptr);
}

template <int dim>
void MGBase<dim>::NonthermalSolves (std::unique_ptr<EquationBase<dim>> equ_ptr)
{}

template <int dim>
void MGBase<dim>::ThermalIterations (std::unique_ptr<EquationBase<dim>> equ_ptr)
{}

template class MGBase<1>
template class MGBase<2>;
template class MGBase<3>;
