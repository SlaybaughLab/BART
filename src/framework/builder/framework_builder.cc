#include "framework/builder/framework_builder.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <sstream>
#include <fstream>

// Convergence classes
#include "convergence/final_checker_or_n.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/moments/multi_moment_checker_max.h"
#include "convergence/parameters/single_parameter_checker.h"

// Domain classes
#include "domain/definition.h"
#include "domain/finite_element/finite_element_gaussian.h"
#include "domain/mesh/mesh_cartesian.h"

// Formulation classes
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"
#include "formulation/stamper.h"
#include "formulation/updater/saaf_updater.h"
#include "formulation/updater/diffusion_updater.h"

// Framework class
#include "framework/framework.h"

// KEffective Updater Classes
#include "calculator/cell/total_aggregated_fission_source.h"
#include "calculator/cell/integrated_fission_source.h"
#include "eigenvalue/k_effective/updater_via_fission_source.h"

// Material classes
#include "material/material_protobuf.h"

// Solver classes
#include "solver/group/single_group_solver.h"
#include "solver/linear/gmres.h"

// Iteration classes
#include "iteration/initializer/initialize_fixed_terms_once.h"
#include "iteration/group/group_solve_iteration.h"
#include "iteration/group/group_source_iteration.h"
#include "iteration/outer/outer_power_iteration.h"
#include "iteration/outer/outer_fixed_source_iteration.h"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

// Results class
#include "results/output_dealii_vtu.h"

// System classes
#include "system/system.h"
#include "system/solution/mpi_group_angular_solution.h"
#include "system/system_functions.h"
#include "system/solution/solution_types.h"

// Instrumentation
#include "instrumentation/builder/instrument_builder.hpp"

namespace bart {

namespace framework {

namespace builder {

namespace  {

using StringColorPair = std::pair<std::string, utility::Color>;

} // namespace

template<int dim>
auto FrameworkBuilder<dim>::BuildFramework(std::string name,
                                           ParametersType& prm)
-> std::unique_ptr<FrameworkType> {

  validator_.Parse(prm);
  // Framework parameters
  int n_angles = 1; // Set to default value of 1 for scalar solve
  const int n_groups = prm.NEnergyGroups();
  const bool need_angular_solution_storage =
      validator_.NeededParts().count(FrameworkPart::AngularSolutionStorage);
  const auto reflective_boundaries = prm.ReflectiveBoundary();
  const bool has_reflective = std::any_of(
      reflective_boundaries.begin(),
      reflective_boundaries.end(),
      [](std::pair<problem::Boundary, bool> pair){ return pair.second; });
  filename_ = prm.OutputFilenameBase();

  instrumentation::builder::InstrumentBuilder instrument_builder;
  using InstrumentName = instrumentation::builder::InstrumentName;

  status_instrument_ptr_ = Shared(
      instrument_builder.BuildInstrument<ColorStatusPair>(
          InstrumentName::kColorStatusToConditionalOstream));

  data_port::StatusDataPort::AddInstrument(status_instrument_ptr_);
  instrumentation::GetPort<data_port::ValidatorStatusPort>(validator_)
      .AddInstrument(status_instrument_ptr_);

  Report("Building framework: " + name + "\n", utility::Color::kGreen);


  auto finite_element_ptr = Shared(BuildFiniteElement(prm));
  auto cross_sections_ptr = Shared(BuildCrossSections(prm));

  auto domain_ptr = Shared(BuildDomain(prm, finite_element_ptr,
                                       ReadMappingFile(prm.MaterialMapFilename())));
  Report("Setting up domain...\n", utility::Color::kReset);
  domain_ptr->SetUpMesh(prm.UniformRefinements()).SetUpDOF();

  // Various objects to be initialized
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr = nullptr;
  UpdaterPointers updater_pointers;
  std::unique_ptr<MomentCalculatorType> moment_calculator_ptr = nullptr;
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solutions_;

  if (prm.TransportModel() != problem::EquationType::kDiffusion) {
    quadrature_set_ptr = BuildQuadratureSet(prm);
    n_angles = quadrature_set_ptr->size();
  };

  if (need_angular_solution_storage) {
    system::SetUpEnergyGroupToAngularSolutionPtrMap(angular_solutions_,
                                                    n_groups, n_angles);
  }


  if (prm.TransportModel() == problem::EquationType::kSelfAdjointAngularFlux) {
    auto stamper_ptr = BuildStamper(domain_ptr);

    auto saaf_formulation_ptr = BuildSAAFFormulation(finite_element_ptr,
                                                     cross_sections_ptr,
                                                     quadrature_set_ptr);
    saaf_formulation_ptr->Initialize(domain_ptr->Cells().at(0));

    if (!has_reflective) {
      updater_pointers = BuildUpdaterPointers(
          std::move(saaf_formulation_ptr),
          std::move(stamper_ptr),
          quadrature_set_ptr);
    } else {
      updater_pointers = BuildUpdaterPointers(
          std::move(saaf_formulation_ptr),
          std::move(stamper_ptr),
          quadrature_set_ptr,
          prm.ReflectiveBoundary(),
          angular_solutions_);
    }
    moment_calculator_ptr = std::move(BuildMomentCalculator(quadrature_set_ptr));

  } else if (prm.TransportModel() == problem::EquationType::kDiffusion) {
    auto diffusion_formulation_ptr = BuildDiffusionFormulation(
        finite_element_ptr,
        cross_sections_ptr);
    diffusion_formulation_ptr->Precalculate(domain_ptr->Cells().at(0));
    auto stamper_ptr = BuildStamper(domain_ptr);

    updater_pointers = BuildUpdaterPointers(
        std::move(diffusion_formulation_ptr),
        std::move(stamper_ptr),
        prm.ReflectiveBoundary());

    moment_calculator_ptr = std::move(BuildMomentCalculator());
  }

  auto initializer_ptr = BuildInitializer(
      updater_pointers.fixed_updater_ptr, n_groups, n_angles);
  auto group_solution_ptr = Shared(BuildGroupSolution(n_angles));
  system::SetUpMPIAngularSolution(*group_solution_ptr, *domain_ptr);

  auto iterative_group_solver_ptr = BuildGroupSolveIteration(
      BuildSingleGroupSolver(),
      BuildMomentConvergenceChecker(1e-6, 10000),
      std::move(moment_calculator_ptr),
      group_solution_ptr,
      updater_pointers,
      BuildMomentMapConvergenceChecker(1e-6, 1000));

  if (need_angular_solution_storage) {
    dynamic_cast<iteration::group::GroupSolveIteration<dim>*>(
        iterative_group_solver_ptr.get())->UpdateThisAngularSolutionMap(
            angular_solutions_);
    validator_.AddPart(FrameworkPart::AngularSolutionStorage);
  };

  std::unique_ptr<OuterIterationType> outer_iteration_ptr;

  if (prm.IsEigenvalueProblem()) {
    outer_iteration_ptr = BuildOuterIteration(
        std::move(iterative_group_solver_ptr),
        BuildParameterConvergenceChecker(1e-6, 10000),
        BuildKEffectiveUpdater(finite_element_ptr, cross_sections_ptr, domain_ptr),
        updater_pointers.fission_source_updater_ptr);
  } else {
    outer_iteration_ptr = BuildOuterIteration(
        std::move(iterative_group_solver_ptr),
        BuildParameterConvergenceChecker(1e-6, 10000));
  };


  auto system_ptr = BuildSystem(n_groups, n_angles, *domain_ptr,
                                group_solution_ptr->solutions().at(0).size(),
                                prm.IsEigenvalueProblem(),
                                need_angular_solution_storage);

  auto results_output_ptr =
      std::make_unique<results::OutputDealiiVtu<dim>>(domain_ptr);

  Validate();

  return std::make_unique<framework::Framework>(
      std::move(system_ptr),
      std::move(initializer_ptr),
      std::move(outer_iteration_ptr),
      std::move(results_output_ptr));
}

// =============================================================================

template<int dim>
auto FrameworkBuilder<dim>::BuildCrossSections(
    const problem::ParametersI& problem_parameters)
    -> std::unique_ptr<CrossSectionType> {
  ReportBuildingComponant("Cross-sections");
  std::unique_ptr<CrossSectionType> return_ptr = nullptr;
  // Default implementation using protocol buffers
  try {
    MaterialProtobuf materials(problem_parameters.MaterialFilenames(),
                               problem_parameters.IsEigenvalueProblem(),
                               problem_parameters.DoNDA(),
                               problem_parameters.NEnergyGroups(),
                               problem_parameters.NumberOfMaterials());
    return_ptr = std::move(std::make_unique<CrossSectionType>(materials));
    ReportBuildSuccess("(default) Cross-sections using protobuf");
  } catch (...) {
    ReportBuildError("(default) Cross-sections using protobuf");
    throw;
  }
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDiffusionFormulation(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const formulation::DiffusionFormulationImpl implementation)
-> std::unique_ptr<DiffusionFormulationType> {
  ReportBuildingComponant("Diffusion formulation");
  std::unique_ptr<DiffusionFormulationType> return_ptr = nullptr;

  if (implementation == formulation::DiffusionFormulationImpl::kDefault) {
    using ReturnType = formulation::scalar::Diffusion<dim>;
    return_ptr = std::move(std::make_unique<ReturnType>(
        finite_element_ptr, cross_sections_ptr));
  }
  ReportBuildSuccess(return_ptr->description());

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDomain(
    ParametersType problem_parameters,
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    std::string material_mapping)
-> std::unique_ptr<DomainType>{
  std::unique_ptr<DomainType> return_ptr = nullptr;

  ReportBuildingComponant("Mesh");
  auto mesh_ptr = std::make_unique<domain::mesh::MeshCartesian<dim>>(
      problem_parameters.SpatialMax(),
      problem_parameters.NCells(),
      material_mapping);
  ReportBuildSuccess(mesh_ptr->description());

  ReportBuildingComponant("Domain");
  return_ptr = std::move(std::make_unique<domain::Definition<dim>>(
      std::move(mesh_ptr),
      finite_element_ptr));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFiniteElement(ParametersType problem_parameters)
-> std::unique_ptr<FiniteElementType>{
  std::unique_ptr<FiniteElementType> return_ptr = nullptr;

  using FiniteElementGaussianType = domain::finite_element::FiniteElementGaussian<dim>;

  ReportBuildingComponant("Cell finite element basis");

  try {
    return_ptr = std::move(std::make_unique<FiniteElementGaussianType>(
        problem::DiscretizationType::kContinuousFEM,
        problem_parameters.FEPolynomialDegree()));

    ReportBuildSuccess(return_ptr->description());
  } catch (...) {
    ReportBuildError();
    throw;
  }
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildUpdaterPointers(
    std::unique_ptr<DiffusionFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::map<problem::Boundary, bool>& reflective_boundaries)
-> UpdaterPointers {
  ReportBuildingComponant("Building Diffusion Formulation updater");
  UpdaterPointers return_struct;

  std::unordered_set<problem::Boundary> reflective_boundary_set;

  for (const auto boundary_pair : reflective_boundaries) {
    if (boundary_pair.second)
      reflective_boundary_set.insert(boundary_pair.first);
  }

  using ReturnType = formulation::updater::DiffusionUpdater<dim>;

  auto diffusion_updater_ptr = std::make_shared<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr),
      reflective_boundary_set);
  ReportBuildSuccess(diffusion_updater_ptr->description());
  return_struct.fixed_updater_ptr = diffusion_updater_ptr;
  return_struct.scattering_source_updater_ptr = diffusion_updater_ptr;
  return_struct.fission_source_updater_ptr = diffusion_updater_ptr;

  return return_struct;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildUpdaterPointers(
    std::unique_ptr<SAAFFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr)
-> UpdaterPointers {
  ReportBuildingComponant("Building SAAF Formulation updater");
  UpdaterPointers return_struct;

  using ReturnType = formulation::updater::SAAFUpdater<dim>;
  auto saaf_updater_ptr = std::make_shared<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr),
      quadrature_set_ptr);
  ReportBuildSuccess(saaf_updater_ptr->description());
  return_struct.fixed_updater_ptr = saaf_updater_ptr;
  return_struct.scattering_source_updater_ptr = saaf_updater_ptr;
  return_struct.fission_source_updater_ptr = saaf_updater_ptr;

  return return_struct;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildUpdaterPointers(
    std::unique_ptr<SAAFFormulationType> formulation_ptr,
    std::unique_ptr<StamperType> stamper_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr,
    const std::map<problem::Boundary, bool>& reflective_boundaries,
    const AngularFluxStorage& angular_flux_storage)
-> UpdaterPointers {
  ReportBuildingComponant("Building SAAF Formulation updater "
                          "(with boundary conditions update)");
  UpdaterPointers return_struct;

  // Transform map into set

  std::unordered_set<problem::Boundary> reflective_boundary_set;

  for (const auto boundary_pair : reflective_boundaries) {
    if (boundary_pair.second)
      reflective_boundary_set.insert(boundary_pair.first);
  }

  using ReturnType = formulation::updater::SAAFUpdater<dim>;
  auto saaf_updater_ptr = std::make_shared<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr),
      quadrature_set_ptr,
      angular_flux_storage,
      reflective_boundary_set);
  ReportBuildSuccess(saaf_updater_ptr->description());

  return_struct.fixed_updater_ptr = saaf_updater_ptr;
  return_struct.scattering_source_updater_ptr = saaf_updater_ptr;
  return_struct.fission_source_updater_ptr = saaf_updater_ptr;
  return_struct.boundary_conditions_updater_ptr = saaf_updater_ptr;

  return return_struct;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildGroupSolveIteration(
    std::unique_ptr<SingleGroupSolverType> single_group_solver_ptr,
    std::unique_ptr<MomentConvergenceCheckerType> moment_convergence_checker_ptr,
    std::unique_ptr<MomentCalculatorType> moment_calculator_ptr,
    const std::shared_ptr<GroupSolutionType>& group_solution_ptr,
    const UpdaterPointers& updater_ptrs,
    std::unique_ptr<MomentMapConvergenceCheckerType> moment_map_convergence_checker_ptr)
    -> std::unique_ptr<GroupSolveIterationType> {
  std::unique_ptr<GroupSolveIterationType> return_ptr = nullptr;

  ReportBuildingComponant("Iterative group solver");

  if (updater_ptrs.boundary_conditions_updater_ptr == nullptr) {
    return_ptr = std::move(
        std::make_unique<iteration::group::GroupSourceIteration<dim>>(
            std::move(single_group_solver_ptr),
            std::move(moment_convergence_checker_ptr),
            std::move(moment_calculator_ptr),
            group_solution_ptr,
            updater_ptrs.scattering_source_updater_ptr,
            std::move(moment_map_convergence_checker_ptr))
    );
  } else {
    return_ptr = std::move(
        std::make_unique<iteration::group::GroupSourceIteration<dim>>(
            std::move(single_group_solver_ptr),
            std::move(moment_convergence_checker_ptr),
            std::move(moment_calculator_ptr),
            group_solution_ptr,
            updater_ptrs.scattering_source_updater_ptr,
            updater_ptrs.boundary_conditions_updater_ptr,
            std::move(moment_map_convergence_checker_ptr))
    );
  }

  validator_.AddPart(FrameworkPart::ScatteringSourceUpdate);
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildGroupSolution(const int n_angles)
-> std::unique_ptr<GroupSolutionType> {
  std::unique_ptr<GroupSolutionType> return_ptr = nullptr;
  ReportBuildingComponant("Group solution");

  return_ptr = std::move(
      std::make_unique<system::solution::MPIGroupAngularSolution>(n_angles));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildInitializer(
    const std::shared_ptr<formulation::updater::FixedUpdaterI>& updater_ptr,
    const int total_groups,
    const int total_angles) -> std::unique_ptr<InitializerType> {
  ReportBuildingComponant("Initializer");

  std::unique_ptr<InitializerType> return_ptr = nullptr;

  using InitializeOnceType = iteration::initializer::InitializeFixedTermsOnce;

  return_ptr = std::move(std::make_unique<InitializeOnceType>(
      updater_ptr, total_groups, total_angles));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildKEffectiveUpdater(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<CrossSectionType>& cross_sections_ptr,
    const std::shared_ptr<DomainType>& domain_ptr)
-> std::unique_ptr<KEffectiveUpdaterType> {
  using AggregatedFissionSource = calculator::cell::TotalAggregatedFissionSource<dim>;
  using IntegratedFissionSource = calculator::cell::IntegratedFissionSource<dim>;
  using ReturnType = eigenvalue::k_effective::UpdaterViaFissionSource;

  ReportBuildingComponant("K_Effective updater");
  std::unique_ptr<KEffectiveUpdaterType> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<ReturnType>(
          std::make_unique<AggregatedFissionSource>(
              std::make_unique<IntegratedFissionSource>(finite_element_ptr,
                                                        cross_sections_ptr),
              domain_ptr),
          2.0,
          10));


  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildMomentCalculator(
    MomentCalculatorImpl implementation)
-> std::unique_ptr<MomentCalculatorType> {
  return BuildMomentCalculator(nullptr, implementation);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildMomentCalculator(
    std::shared_ptr<QuadratureSetType> quadrature_set_ptr,
    FrameworkBuilder::MomentCalculatorImpl implementation)
-> std::unique_ptr<MomentCalculatorType> {
  ReportBuildingComponant("Moment calculator");
  std::unique_ptr<MomentCalculatorType> return_ptr = nullptr;

  try {
    return_ptr = std::move(quadrature::factory::MakeMomentCalculator<dim>(
        implementation, quadrature_set_ptr));

    if (implementation == MomentCalculatorImpl::kScalarMoment) {
      ReportBuildSuccess("(default) calculator for scalar solve");
    } else if (implementation == MomentCalculatorImpl::kZerothMomentOnly) {
      ReportBuildSuccess("(default) calculator for 0th moment only");
    } else {
      AssertThrow(false,
                  dealii::ExcMessage("Unsupported implementation of moment "
                                     "calculator specified in call to "
                                     "BuildMomentCalculator"))
    }
  } catch (...) {
    ReportBuildError();
    throw;
  }

  return return_ptr;
}



template<int dim>
auto FrameworkBuilder<dim>::BuildMomentConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentConvergenceCheckerType>{
  //TODO(Josh): Add option for using other than L1Norm
  ReportBuildingComponant("Moment convergence checker");

  using CheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  using FinalCheckerType = convergence::FinalCheckerOrN<
      system::moments::MomentVector,
      convergence::moments::SingleMomentCheckerI>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));
  return_ptr->SetMaxIterations(max_iterations);
  return return_ptr;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildMomentMapConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentMapConvergenceCheckerType> {
  ReportBuildingComponant("Moment map convergence checker");

  using SingleCheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  using CheckerType = convergence::moments::MultiMomentCheckerMax;
  using FinalCheckerType = convergence::FinalCheckerOrN<
      const system::moments::MomentsMap,
      convergence::moments::MultiMomentCheckerI>;
  auto checker_ptr = std::make_unique<CheckerType>(
      std::make_unique<SingleCheckerType>(max_delta));
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(checker_ptr));
  return_ptr->SetMaxIterations(max_iterations);
  return return_ptr;
}
template <int dim>
auto FrameworkBuilder<dim>::BuildOuterIteration(
    std::unique_ptr<GroupSolveIterationType> group_iteration_ptr,
    std::unique_ptr<ParameterConvergenceCheckerType> convergence_checker_ptr)
    -> std::unique_ptr<OuterIterationType> {
  ReportBuildingComponant("Outer iteration");
  std::unique_ptr<OuterIterationType> return_ptr = nullptr;
  using ReturnType = iteration::outer::OuterFixedSourceIteration;

  return_ptr = std::move(std::make_unique<ReturnType>(
      std::move(group_iteration_ptr),
      std::move(convergence_checker_ptr)));

  using ConvergenceDataPort = iteration::outer::data_names::ConvergenceStatusPort;
  using StatusPort =  iteration::outer::data_names::StatusPort;

//  instrumentation::GetPort<ConvergenceDataPort>(*return_ptr)
//      .AddInstrument(Shared(BuildConvergenceInstrument()));
//  instrumentation::GetPort<StatusPort>(*return_ptr)
//      .AddInstrument(status_instrument_ptr_);

  ReportBuildSuccess(return_ptr->description());

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildOuterIteration(
    std::unique_ptr<GroupSolveIterationType> group_solve_iteration_ptr,
    std::unique_ptr<ParameterConvergenceCheckerType> parameter_convergence_checker_ptr,
    std::unique_ptr<KEffectiveUpdaterType> k_effective_updater_ptr,
    const std::shared_ptr<FissionSourceUpdaterType>& fission_source_updater_ptr)
-> std::unique_ptr<OuterIterationType> {
  std::unique_ptr<OuterIterationType> return_ptr = nullptr;
  ReportBuildingComponant("Outer Iteration");

  using DefaultOuterPowerIteration = iteration::outer::OuterPowerIteration;

  return_ptr = std::move(
      std::make_unique<DefaultOuterPowerIteration>(
          std::move(group_solve_iteration_ptr),
          std::move(parameter_convergence_checker_ptr),
          std::move(k_effective_updater_ptr),
          fission_source_updater_ptr));

  using ConvergenceDataPort = iteration::outer::data_names::ConvergenceStatusPort;
  using StatusPort =  iteration::outer::data_names::StatusPort;
  using IterationErrorPort = iteration::outer::data_names::IterationErrorPort;

//  instrumentation::GetPort<ConvergenceDataPort>(*return_ptr)
//      .AddInstrument(Shared(BuildConvergenceInstrument()));
//  instrumentation::GetPort<StatusPort>(*return_ptr)
//      .AddInstrument(Shared(BuildStatusInstrument()));
//  instrumentation::GetPort<IterationErrorPort>(*return_ptr)
//      .AddInstrument(Shared(BuildIterationErrorInstrument(
//          filename_ + "_iteration_error.csv")));

  validator_.AddPart(FrameworkPart::FissionSourceUpdate);
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildParameterConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<ParameterConvergenceCheckerType>{
  ReportBuildingComponant("Parameter (double) convergence checker");
  using CheckerType = convergence::parameters::SingleParameterChecker;
  using FinalCheckerType = convergence::FinalCheckerOrN<double, CheckerType>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));
  return_ptr->SetMaxIterations(max_iterations);

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildQuadratureSet(ParametersType problem_parameters)
-> std::shared_ptr<QuadratureSetType> {
  ReportBuildingComponant("quadrature set");
  using QuadratureGeneratorType = quadrature::QuadratureGeneratorI<dim>;

  std::shared_ptr<QuadratureSetType> return_ptr = nullptr;
  std::shared_ptr<QuadratureGeneratorType > quadrature_generator_ptr = nullptr;

  const int order_value = problem_parameters.AngularQuadOrder();
  switch (problem_parameters.AngularQuad()) {
    case problem::AngularQuadType::kLevelSymmetricGaussian: {
      AssertThrow(dim == 3, dealii::ExcMessage("Error in BuildQuadratureSet "
                                               "LSGC only available for 3D"))
      quadrature_generator_ptr =
          quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
              quadrature::Order(order_value),
              quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
      break;
    }
    case problem::AngularQuadType::kGaussLegendre: {
      AssertThrow(dim == 1, dealii::ExcMessage("Error in BuildQuadratureSet "
                                               "GaussLegendre only available "
                                               "for 1D"))
      quadrature_generator_ptr =
          quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
              quadrature::Order(order_value),
              quadrature::AngularQuadratureSetType::kGaussLegendre);
      break;
    }
    default: {
        AssertThrow(false,
                    dealii::ExcMessage("No supported quadratures for this dimension "
                                       "and transport model"))
      break;
    }
  }
  ReportBuildSuccess(quadrature_generator_ptr->description());

  return_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

  auto quadrature_points = quadrature::utility::GenerateAllPositiveX<dim>(
      quadrature_generator_ptr->GenerateSet());

  quadrature::factory::FillQuadratureSet<dim>(return_ptr.get(),
                                              quadrature_points);

  return return_ptr;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildSAAFFormulation(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr,
    const formulation::SAAFFormulationImpl implementation)
-> std::unique_ptr<SAAFFormulationType> {
  ReportBuildingComponant("Building SAAF Formulation");
  std::unique_ptr<SAAFFormulationType> return_ptr;

  if (implementation == formulation::SAAFFormulationImpl::kDefault) {
    using ReturnType = formulation::angular::SelfAdjointAngularFlux<dim>;

    return_ptr = std::move(std::make_unique<ReturnType>(finite_element_ptr,
                                                        cross_sections_ptr,
                                                        quadrature_set_ptr));
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSingleGroupSolver(
    const int max_iterations,
    const double convergence_tolerance)
-> std::unique_ptr<SingleGroupSolverType> {
  ReportBuildingComponant("Single group solver");
  std::unique_ptr<SingleGroupSolverType> return_ptr = nullptr;

  auto linear_solver_ptr = std::make_unique<solver::GMRES>(max_iterations,
                                                           convergence_tolerance);
  ReportBuildSuccess("GMRES: tol = " + std::to_string(convergence_tolerance)
                            + "iter_max = " + std::to_string(max_iterations));
  return_ptr = std::move(std::make_unique<solver::group::SingleGroupSolver>(
          std::move(linear_solver_ptr)));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSystem(
    const int total_groups,
    const int total_angles,
    const DomainType& domain,
    const std::size_t solution_size,
    bool is_eigenvalue_problem,
    bool need_rhs_boundary_condition) -> std::unique_ptr<SystemType> {
  std::unique_ptr<SystemType> return_ptr;

  ReportBuildingComponant("system");
  try {
    return_ptr = std::move(std::make_unique<SystemType>());
    system::InitializeSystem(*return_ptr, total_groups, total_angles,
                             is_eigenvalue_problem, need_rhs_boundary_condition);
    system::SetUpSystemTerms(*return_ptr, domain);
    system::SetUpSystemMoments(*return_ptr, solution_size);
    ReportBuildSuccess("system");
  } catch (...) {
    ReportBuildError("system initialization error.");
    throw;
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildStamper(
    const std::shared_ptr<DomainType>& domain_ptr)
-> std::unique_ptr<StamperType> {
  ReportBuildingComponant("Stamper");
  std::unique_ptr<StamperType> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<formulation::Stamper<dim>>(domain_ptr));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}
template<int dim>
std::string FrameworkBuilder<dim>::ReadMappingFile(std::string filename) {
  ReportBuildingComponant("Reading mapping file: ");

  std::ifstream mapping_file(filename);
  if (mapping_file.is_open()) {
    ReportBuildSuccess(filename);

    return std::string(
        (std::istreambuf_iterator<char>(mapping_file)),
        std::istreambuf_iterator<char>());
  } else {
    ReportBuildError("Error reading " + filename);
    AssertThrow(false,
                dealii::ExcMessage("Failed to open material mapping file"))
  }
}

template<int dim>
void FrameworkBuilder<dim>::Validate() const {
  validator_.ReportValidation();
}

template class FrameworkBuilder<1>;
template class FrameworkBuilder<2>;
template class FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart
