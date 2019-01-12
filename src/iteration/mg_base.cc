#include "mg_base.h"

#include "gauss_seidel.h"

template <int dim>
std::unique_ptr<MGBase<dim>> MGBase<dim>::CreateMGIteration (
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr) {

  const std::string mg_iteration_name(prm.get("mg solver name"));

  std::unordered_map<std::string, MGIterationType> mg_iteration_name_map =
      {{"gs", MGIterationType::kGaussSeidel}};

  std::unique_ptr<MGBase<dim>> iteration_ptr;

  switch(mg_iteration_name_map[mg_iteration_name]) {
    case MGIterationType::kGaussSeidel: {
      iteration_ptr.reset(new GaussSeidel<dim>(prm, dat_ptr));
      break;
    }
    default: {
      assert(false);
      break;
    }
  }
  return std::move(iteration_ptr);
}

template <int dim>
MGBase<dim>::MGBase (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    IterationBase<dim> (prm, dat_ptr),
    g_thermal_(prm.get_integer("thermal group boundary")),
    ig_ptr_(IGBase<dim>::CreateIGIteration(prm, dat_ptr)),
    err_phi_tol_(1.0e-5) {
  AssertThrow (g_thermal_<this->n_group_ && g_thermal_>=0,
               dealii::ExcMessage("Invalid thermal upper boundary"));
}

template <int dim>
void MGBase<dim>::DoIterations (std::unordered_map<std::string,
    std::unique_ptr<EquationBase<dim>>> &equ_ptrs) {
  // this function will be called only when it's not an eigenvalue problem. Copy
  // the following assertion in ALL derived class of MGBase.
  AssertThrow (!this->is_eigen_problem_,
               dealii::ExcMessage("This function shall not be called if it's eigen prob."));
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
void MGBase<dim>::MGIterations (std::unique_ptr<EquationBase<dim>> &equ_ptr) {
  // this function needs to be overridden if JFNK is desired
  // solve for epithermal and fast groups
  NonthermalSolves (equ_ptr);
  // thermal group iterations
  ThermalIterations (equ_ptr);
}

template <int dim>
void MGBase<dim>::NonthermalSolves (std::unique_ptr<EquationBase<dim>> &)
{}

template <int dim>
void MGBase<dim>::ThermalIterations (std::unique_ptr<EquationBase<dim>> &)
{}

template class MGBase<1>;
template class MGBase<2>;
template class MGBase<3>;
