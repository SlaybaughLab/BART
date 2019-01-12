#include "ig_base.h"

#include "source_iteration.h"

template <int dim>
IGBase<dim>::IGBase (const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr)
    :
    IterationBase<dim> (prm, dat_ptr),
    err_phi_tol_(1.0e-6) {}

template <int dim>
std::unique_ptr<IGBase<dim>> IGBase<dim>::CreateIGIteration (
      const dealii::ParameterHandler &prm,
      std::shared_ptr<FundamentalData<dim>> &dat_ptr) {
  const std::string ig_iteration_name(prm.get("in group solver name"));
  std::unordered_map<std::string, IGIterationType> ig_iteration_name_map =
      {{"si", IGIterationType::SourceIteration}};

  std::unique_ptr<IGBase<dim>> iteration_ptr;

  switch(ig_iteration_name_map[ig_iteration_name]) {
    case IGIterationType::SourceIteration: {
      iteration_ptr.reset(new SourceIteration<dim>(prm, dat_ptr));
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
void IGBase<dim>::IGIterations (std::unique_ptr<EquationBase<dim>> &equ_ptr,
    const int &g) {
  // the default is for diffusion like system, SPN and PN
  equ_ptr->SolveInGroup (g);
}

template class IGBase<1>;
template class IGBase<2>;
template class IGBase<3>;
