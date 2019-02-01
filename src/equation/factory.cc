#include "factory.h"

#include "../problem/parameter_types.h"
#include "../problem/locator.h"
#include "diffusion.h"
#include "self_adjoint_angular_flux.h"
#include "even_parity.h"

namespace bart {

namespace equation {

template <int dim>
std::unique_ptr<EquationBase<dim>> Factory<dim>::CreateEquation(
    const dealii::ParameterHandler &prm,
    std::shared_ptr<FundamentalData<dim>> &dat_ptr) {

  using EquationType = bart::problem::EquationType;
  
  // Get transport model
  auto problem_parameters = bart::problem::Locator::GetParameters();
  const std::string equation_name(prm.get("transport model"));
  
  const EquationType equation_type{problem_parameters->TransportModel()};

  std::unique_ptr<EquationBase<dim>> eq_ptr;

  switch(equation_type) {
    case EquationType::kSelfAdjointAngularFlux: {
      eq_ptr.reset(
          new SelfAdjointAngularFlux<dim>(equation_name, prm, dat_ptr));
      break;
    }
    case EquationType::kEvenParity: {
      eq_ptr.reset(new EvenParity<dim>(equation_name, prm, dat_ptr));
      break;
    }
    case EquationType::kDiffusion: {
      eq_ptr.reset(new Diffusion<dim>(equation_name, prm, dat_ptr));
      break;
    }
    default: {
      assert(false);
    }      
  }
  
  return std::move(eq_ptr);
}

template <int dim>
std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>>
GetEquations(const dealii::ParameterHandler &prm,
             std::shared_ptr<FundamentalData<dim>> &dat_ptr) {

  std::unordered_map<std::string, std::unique_ptr<EquationBase<dim>>>
      return_map;

  return_map[prm.get("transport model")] =
      EquationBase<dim>::CreateEquation(prm, dat_ptr);
  
  return std::move(return_map);
}

template 
std::unordered_map<std::string, std::unique_ptr<EquationBase<1>>>
GetEquations(const dealii::ParameterHandler &prm,
             std::shared_ptr<FundamentalData<1>> &dat_ptr);

template 
std::unordered_map<std::string, std::unique_ptr<EquationBase<2>>>
GetEquations(const dealii::ParameterHandler &prm,
             std::shared_ptr<FundamentalData<2>> &dat_ptr);

template 
std::unordered_map<std::string, std::unique_ptr<EquationBase<3>>>
GetEquations(const dealii::ParameterHandler &prm,
             std::shared_ptr<FundamentalData<3>> &dat_ptr);


} // namespace equation

} // namespace bart
