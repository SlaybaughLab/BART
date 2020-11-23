#include "framework/framework_helper.hpp"

#include "framework/builder/framework_validator.hpp"
#include "system/system_helper.hpp"
#include "system/solution/solution_types.h"

#include <fmt/color.h>
#include <system/system_helper.hpp>

namespace bart::framework {

namespace  {

template <typename T> inline std::shared_ptr<T> Shared(std::unique_ptr<T> to_convert_ptr) { return to_convert_ptr; }
} // namespace

template<int dim>
FrameworkHelper<dim>::FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr)
    : system_helper_ptr_(system_helper_ptr) {
  AssertThrow(system_helper_ptr_ != nullptr, dealii::ExcMessage("Error in constructor of framework helper, system"
                                                                "helper pointer passed is null"))
}

template<int dim>
auto FrameworkHelper<dim>::BuildFramework(
    builder::FrameworkBuilderI<dim>& builder,
    const framework::FrameworkParameters& parameters) -> std::unique_ptr<framework::FrameworkI> {
  using FrameworkPart = framework::builder::FrameworkPart;
  using MomentCalculator = typename builder::FrameworkBuilderI<dim>::MomentCalculator;
  using MomentCalculatorImpl = typename builder::FrameworkBuilderI<dim>::MomentCalculatorImpl;
  using UpdaterPointers = typename builder::FrameworkBuilderI<dim>::UpdaterPointers;
  using QuadratureSet = typename builder::FrameworkBuilderI<dim>::QuadratureSet;
  using Validator = framework::builder::FrameworkValidator;

  Validator validator;
  validator.Parse(parameters);

  const int n_groups{ parameters.neutron_energy_groups };
  const bool need_angular_solution_storage{ validator.NeededParts().contains(FrameworkPart::AngularSolutionStorage) };
  const bool has_reflective_boundaries { !parameters.reflective_boundaries.empty() };
  const std::string output_filename_base { parameters.output_filename_base };

  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "Building framework: {}\n", parameters.name);

  auto finite_element_ptr = Shared(builder.BuildFiniteElement(parameters.cell_finite_element_type,
                                                              parameters.discretization_type,
                                                              parameters.polynomial_degree));
  auto domain_ptr = Shared(builder.BuildDomain(parameters.domain_size, parameters.number_of_cells,
                                               finite_element_ptr, parameters.material_mapping));

  fmt::print("Setting up domain...\n");
  domain_ptr->SetUpMesh(parameters.uniform_refinements);
  domain_ptr->SetUpDOF();

  // These objects will be set up differently depending on the implementation
  std::shared_ptr<QuadratureSet> quadrature_set_ptr{ nullptr };
  UpdaterPointers updater_pointers;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr{ nullptr };


  // Set up for Angular/Scalar solve
  if (parameters.equation_type != problem::EquationType::kDiffusion) {
    // Angular solve
    AssertThrow(parameters.angular_quadrature_order.has_value(),
                dealii::ExcMessage("Error building framework, equation type requires quadrature but order is null"))
    quadrature_set_ptr = builder.BuildQuadratureSet(parameters.angular_quadrature_type,
                                                    parameters.angular_quadrature_order.value());
    moment_calculator_ptr = builder.BuildMomentCalculator(quadrature_set_ptr, MomentCalculatorImpl::kZerothMomentOnly);
  } else {
    // Scalar solve
    moment_calculator_ptr = builder.BuildMomentCalculator(MomentCalculatorImpl::kScalarMoment);
  }

  const int n_angles { quadrature_set_ptr == nullptr ? 1 : static_cast<int>(quadrature_set_ptr->size()) };

  // Set up angular solutions if needed
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solutions_;
  if (need_angular_solution_storage)
    system_helper_ptr_->SetUpEnergyGroupToAngularSolutionPtrMap(angular_solutions_, n_groups, n_angles);

  //TODO: Add overload that makes this unecessary
  std::map<problem::Boundary, bool> reflective_boundaries {
      {problem::Boundary::kXMin, false}, {problem::Boundary::kXMax, false},
      {problem::Boundary::kYMin, false}, {problem::Boundary::kYMax, false},
      {problem::Boundary::kZMin, false}, {problem::Boundary::kZMax, false},
  };
  for (auto& boundary : parameters.reflective_boundaries) {
    reflective_boundaries.at(boundary) = true;
  }

  // Formulation specific builds
  if (parameters.equation_type == problem::EquationType::kSelfAdjointAngularFlux) {
    auto saaf_formulation_ptr = builder.BuildSAAFFormulation(finite_element_ptr,
                                                             parameters.cross_sections_.value(),
                                                             quadrature_set_ptr,
                                                             formulation::SAAFFormulationImpl::kDefault);
    saaf_formulation_ptr->Initialize(domain_ptr->Cells().at(0));
    if (has_reflective_boundaries) {
      updater_pointers = builder.BuildUpdaterPointers(std::move(saaf_formulation_ptr),
                                                      builder.BuildStamper(domain_ptr),
                                                      quadrature_set_ptr,
                                                      reflective_boundaries,
                                                      angular_solutions_);
    } else {
      updater_pointers = builder.BuildUpdaterPointers(std::move(saaf_formulation_ptr),
                                                      builder.BuildStamper(domain_ptr),
                                                      quadrature_set_ptr);
    }
  } else if (parameters.equation_type == problem::EquationType::kDiffusion) {
    auto diffusion_formulation_ptr = builder.BuildDiffusionFormulation(finite_element_ptr,
                                                                       parameters.cross_sections_.value(),
                                                                       formulation::DiffusionFormulationImpl::kDefault);
    diffusion_formulation_ptr->Precalculate(domain_ptr->Cells().at(0));
    updater_pointers = builder.BuildUpdaterPointers(std::move(diffusion_formulation_ptr),
                                                    builder.BuildStamper(domain_ptr),
                                                    reflective_boundaries);
  }

  auto initializer_ptr = builder.BuildInitializer(updater_pointers.fixed_updater_ptr,
                                                  parameters.neutron_energy_groups,
                                                  n_angles);

  auto group_solution_ptr = Shared(builder.BuildGroupSolution(n_angles));
  system_helper_ptr_->SetUpMPIAngularSolution(*group_solution_ptr, *domain_ptr, 1.0);

  auto group_iteration_ptr = builder.BuildGroupSolveIteration(
      builder.BuildSingleGroupSolver(10000, 1e-10),
      builder.BuildMomentConvergenceChecker(1e-6, 10000),
      std::move(moment_calculator_ptr),
      group_solution_ptr,
      updater_pointers,
      builder.BuildMomentMapConvergenceChecker(1e-6, 1000));

  if (need_angular_solution_storage) {
    group_iteration_ptr->UpdateThisAngularSolutionMap(angular_solutions_);
    validator.AddPart(FrameworkPart::AngularSolutionStorage);
  }

  if (parameters.eigen_solver_type.has_value()){
    auto k_effective_updater = builder.BuildKEffectiveUpdater(finite_element_ptr,
                                                              parameters.cross_sections_.value(),
                                                              domain_ptr);
    auto parameter_convergence_checker_ptr = builder.BuildParameterConvergenceChecker(1e-6, 1000);
  } else {
    auto parameter_convergence_checker_ptr = builder.BuildParameterConvergenceChecker(1e-6, 1000);
  }

  return nullptr;
}

template class FrameworkHelper<1>;
template class FrameworkHelper<2>;
template class FrameworkHelper<3>;

} // namespace bart::framework
