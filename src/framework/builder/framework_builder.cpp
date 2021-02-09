#include "framework/builder/framework_builder.hpp"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <sstream>
#include <fstream>
#include <quadrature/calculators/angular_flux_integrator.hpp>
#include <quadrature/calculators/quadrature_calculators_factories.hpp>
#include <calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp>
#include <calculator/drift_diffusion/factory.hpp>
#include <formulation/scalar/scalar_formulation_factory.hpp>
#include <formulation/updater/formulation_updater_factories.hpp>
#include <iteration/subroutine/get_scalar_flux_from_framework.hpp>

// Builders & factories
#include "solver/builder/solver_builder.hpp"

// Convergence classes
#include "convergence/final_checker_or_n.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/moments/multi_moment_checker_max.h"
#include "convergence/parameters/single_parameter_checker.h"

// Domain classes
#include "domain/definition.h"
#include "domain/finite_element/finite_element_gaussian.hpp"
#include "domain/mesh/mesh_cartesian.h"

// Formulation classes
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"
#include "formulation/scalar/drift_diffusion.hpp"
#include "formulation/stamper.h"
#include "formulation/updater/saaf_updater.h"
#include "formulation/updater/diffusion_updater.hpp"
#include "formulation/updater/drift_diffusion_updater.hpp"

// Framework class
#include "framework/framework.hpp"

// KEffective Updater Classes
#include "calculator/cell/total_aggregated_fission_source.hpp"
#include "calculator/cell/integrated_fission_source.hpp"
#include "eigenvalue/k_effective/updater_via_fission_source.h"
#include "eigenvalue/k_effective/updater_via_rayleigh_quotient.hpp"

// Material classes
#include "material/material_protobuf.h"

// Iteration classes
#include "iteration/initializer/initialize_fixed_terms_once.h"
#include "iteration/initializer/initialize_fixed_terms_reset_moments.hpp"
#include "iteration/group/group_solve_iteration.h"
#include "iteration/group/group_source_iteration.h"
#include "iteration/outer/outer_power_iteration.hpp"
#include "iteration/outer/outer_fixed_source_iteration.hpp"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

// Results class
#include "results/output_dealii_vtu.h"

// System classes
#include "system/system.hpp"
#include "system/solution/mpi_group_angular_solution.h"
#include "system/solution/solution_types.h"

// Instrumentation
#include "instrumentation/builder/instrument_builder.hpp"

namespace bart::framework::builder {

namespace  {

using InstrumentBuilder = instrumentation::builder::InstrumentBuilder;
using InstrumentName = instrumentation::builder::InstrumentName;
using StringColorPair = std::pair<std::string, utility::Color>;

} // namespace

template<int dim>
FrameworkBuilder<dim>::FrameworkBuilder(std::unique_ptr<Validator> validator_ptr)
    : validator_ptr_(std::move(validator_ptr)) {}

// =============================================================================

template<int dim>
auto FrameworkBuilder<dim>::BuildAngularFluxIntegrator(const std::shared_ptr<QuadratureSet> quadrature_set_ptr)
-> std::unique_ptr<AngularFluxIntegrator> {
  return quadrature::calculators::AngularFluxIntegrator<dim>::Factory::get()
      .GetConstructor(quadrature::calculators::AngularFluxIntegratorName::kDefaultImplementation)(quadrature_set_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDiffusionFormulation(const std::shared_ptr<FiniteElement>& finite_element_ptr,
                                                      const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
                                                      const DiffusionFormulationImpl implementation)
-> std::unique_ptr<DiffusionFormulation> {
  ReportBuildingComponant("Diffusion formulation");
  std::unique_ptr<DiffusionFormulation> return_ptr = nullptr;

  if (implementation == DiffusionFormulationImpl::kDefault) {
    using ReturnType = formulation::scalar::Diffusion<dim>;
    return_ptr = std::move(std::make_unique<ReturnType>(finite_element_ptr, cross_sections_ptr));
  }
  ReportBuildSuccess(return_ptr->description());

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDriftDiffusionFormulation(
    const std::shared_ptr<AngularFluxIntegrator>& angular_flux_integrator_ptr,
    const std::shared_ptr<FiniteElement>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr) -> std::unique_ptr<DriftDiffusionFormulation> {
  auto drift_diffusion_vector_calculator_ptr = Shared(calculator::drift_diffusion::DriftDiffusionVectorCalculatorIFactory<dim>::get()
      .GetConstructor(calculator::drift_diffusion::DriftDiffusionVectorCalculatorName::kDefaultImplementation)());
  return formulation::scalar::DriftDiffusion<dim>::Factory::get()
      .GetConstructor(formulation::scalar::DriftDiffusionFormulationName::kDefaultImplementation)(
          finite_element_ptr, cross_sections_ptr, drift_diffusion_vector_calculator_ptr, angular_flux_integrator_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildDomain(FrameworkParameters::DomainSize domain_size,
                                        FrameworkParameters::NumberOfCells number_of_cells,
                                        const std::shared_ptr<FiniteElement>& finite_element_ptr,
                                        std::string material_mapping) -> std::unique_ptr<Domain> {
  std::unique_ptr<Domain> return_ptr = nullptr;
  try {
    ReportBuildingComponant("Mesh");
    auto mesh_ptr = std::make_unique<domain::mesh::MeshCartesian<dim>>(
        domain_size.get(), number_of_cells.get(), material_mapping);
    ReportBuildSuccess(mesh_ptr->description());

    ReportBuildingComponant("Domain");
    return_ptr = std::move(std::make_unique<domain::Definition<dim>>(std::move(mesh_ptr), finite_element_ptr));
    ReportBuildSuccess(return_ptr->description());
  } catch (...) {
    ReportBuildError();
    throw;
  }
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildFiniteElement(problem::CellFiniteElementType finite_element_type,
                                               problem::DiscretizationType discretization_type,
                                               FrameworkParameters::PolynomialDegree polynomial_degree)
-> std::unique_ptr<FiniteElement> {
  ReportBuildingComponant("Cell finite element basis");
  std::unique_ptr<FiniteElement> return_ptr{ nullptr };

  try {
    AssertThrow(polynomial_degree.get() > 0, dealii::ExcMessage("Bad polynomial degree"))
    switch (finite_element_type) {
      case problem::CellFiniteElementType::kGaussian: {
        using ReturnType = domain::finite_element::FiniteElementGaussian<dim>;
        return_ptr = std::move(std::make_unique<ReturnType>(discretization_type, polynomial_degree.get()));
      }
    }
  } catch (...) {
    ReportBuildError();
    throw;
  }
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildUpdaterPointers(
    std::unique_ptr<DiffusionFormulation> diffusion_formulation_ptr,
    std::unique_ptr<DriftDiffusionFormulation> drift_diffusion_formulation_ptr,
    std::shared_ptr<Stamper> stamper_ptr,
    std::shared_ptr<AngularFluxIntegrator> angular_flux_integrator_ptr,
    std::shared_ptr<SphericalHarmonicMoments> higher_order_moments_ptr,
    AngularFluxStorage& angular_flux_storage,
    const std::map<problem::Boundary, bool>& reflective_boundaries) -> UpdaterPointers {
  ReportBuildingComponant("Building Drift-Diffusion Formulation updater");
  UpdaterPointers return_struct;

  std::unordered_set<problem::Boundary> reflective_boundary_set;

  for (const auto boundary_pair : reflective_boundaries) {
    if (boundary_pair.second)
      reflective_boundary_set.insert(boundary_pair.first);
  }

  using ReturnType = formulation::updater::DriftDiffusionUpdater<dim>;
  auto implementation_name{ formulation::updater::DriftDiffusionUpdaterName::kDefaultImplementation };
  auto drift_diffusion_updater_ptr = Shared(ReturnType::Factory::get().GetConstructor(implementation_name)(
      std::move(diffusion_formulation_ptr),
      std::move(drift_diffusion_formulation_ptr),
      stamper_ptr,
      angular_flux_integrator_ptr,
      higher_order_moments_ptr,
      angular_flux_storage,
      reflective_boundary_set));

  return_struct.fixed_updater_ptr = drift_diffusion_updater_ptr;
  return_struct.fission_source_updater_ptr = drift_diffusion_updater_ptr;
  return_struct.scattering_source_updater_ptr = drift_diffusion_updater_ptr;

  return return_struct;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildUpdaterPointers(
    std::unique_ptr<DiffusionFormulation> formulation_ptr,
    std::unique_ptr<Stamper> stamper_ptr,
    const std::map<problem::Boundary, bool>& reflective_boundaries) -> UpdaterPointers {
  ReportBuildingComponant("Building Diffusion Formulation updater");
  UpdaterPointers return_struct;

  std::unordered_set<problem::Boundary> reflective_boundary_set;

  for (const auto boundary_pair : reflective_boundaries) {
    if (boundary_pair.second)
      reflective_boundary_set.insert(boundary_pair.first);
  }

  using ReturnType = formulation::updater::DiffusionUpdater<dim>;

  auto diffusion_updater_ptr = std::make_shared<ReturnType>(std::move(formulation_ptr),
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
    std::unique_ptr<SAAFFormulation> formulation_ptr,
    std::unique_ptr<Stamper> stamper_ptr,
    const std::shared_ptr<QuadratureSet>& quadrature_set_ptr) -> UpdaterPointers {
  ReportBuildingComponant("Building SAAF Formulation updater");
  UpdaterPointers return_struct;

  using ReturnType = formulation::updater::SAAFUpdater<dim>;
  auto saaf_updater_ptr = std::make_shared<ReturnType>(std::move(formulation_ptr),
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
    std::unique_ptr<SAAFFormulation> formulation_ptr,
    std::unique_ptr<Stamper> stamper_ptr,
    const std::shared_ptr<QuadratureSet>& quadrature_set_ptr,
    const std::map<problem::Boundary, bool>& reflective_boundaries,
    const AngularFluxStorage& angular_flux_storage) -> UpdaterPointers {
  ReportBuildingComponant("Building SAAF Formulation updater (with boundary conditions update)");
  UpdaterPointers return_struct;

  // Transform map into set

  std::unordered_set<problem::Boundary> reflective_boundary_set;

  for (const auto boundary_pair : reflective_boundaries) {
    if (boundary_pair.second)
      reflective_boundary_set.insert(boundary_pair.first);
  }

  using ReturnType = formulation::updater::SAAFUpdater<dim>;
  auto saaf_updater_ptr = std::make_shared<ReturnType>(std::move(formulation_ptr),
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
    std::unique_ptr<SingleGroupSolver> single_group_solver_ptr,
    std::unique_ptr<MomentConvergenceChecker> moment_convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    const std::shared_ptr<GroupSolution>& group_solution_ptr,
    const UpdaterPointers& updater_ptrs,
    std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr)
    -> std::unique_ptr<GroupSolveIteration> {
  std::unique_ptr<GroupSolveIteration> return_ptr = nullptr;

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

  using ConvergenceStatusPort = iteration::group::data_ports::ConvergenceStatusPort;
  using StatusPort = iteration::group::data_ports::StatusPort;

  instrumentation::GetPort<ConvergenceStatusPort>(*return_ptr)
      .AddInstrument(convergence_status_instrument_ptr_);
  instrumentation::GetPort<StatusPort>(*return_ptr)
      .AddInstrument(status_instrument_ptr_);

  validator_ptr_->AddPart(FrameworkPart::ScatteringSourceUpdate);
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildGroupSolution(const int n_angles) -> std::unique_ptr<GroupSolution> {
  std::unique_ptr<GroupSolution> return_ptr = nullptr;
  ReportBuildingComponant("Group solution");

  return_ptr = std::move(std::make_unique<system::solution::MPIGroupAngularSolution>(n_angles));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildInitializer(const std::shared_ptr<FixedTermUpdater>& updater_ptr,
                                             const int total_groups,
                                             const int total_angles,
                                             const InitializerName implementation) -> std::unique_ptr<Initializer> {
  ReportBuildingComponant("Initializer");

  std::unique_ptr<Initializer> return_ptr = nullptr;

  if (implementation == InitializerName::kInitializeFixedTermsOnce) {
    using InitializeOnceType = iteration::initializer::InitializeFixedTermsOnce;
    return_ptr = std::move(std::make_unique<InitializeOnceType>(updater_ptr, total_groups, total_angles));
  } else if (implementation == InitializerName::kInitializeFixedTermsAndResetMoments) {
    using InitializeAndResetMoments = iteration::initializer::InitializeFixedTermsResetMoments;
    return_ptr = std::move(std::make_unique<InitializeAndResetMoments>(updater_ptr, total_groups, total_angles));
  }

  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildInitializer(const std::shared_ptr<FixedTermUpdater>& updater_ptr,
                                             const int total_groups,
                                             const int total_angles) -> std::unique_ptr<Initializer> {
  return BuildInitializer(updater_ptr, total_groups, total_angles, InitializerName::kInitializeFixedTermsOnce);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildKEffectiveUpdater() -> std::unique_ptr<KEffectiveUpdater> {
  using ReturnType = eigenvalue::k_effective::UpdaterViaRayleighQuotient;
  ReportBuildingComponant("K_Effective updater");
  std::unique_ptr<KEffectiveUpdater> return_ptr{ nullptr };
  return_ptr = std::move(std::make_unique<ReturnType>());
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildKEffectiveUpdater(
    const std::shared_ptr<FiniteElement>& finite_element_ptr,
    const std::shared_ptr<CrossSections>& cross_sections_ptr,
    const std::shared_ptr<Domain>& domain_ptr)
-> std::unique_ptr<KEffectiveUpdater> {
  using AggregatedFissionSource = calculator::cell::TotalAggregatedFissionSource<dim>;
  using IntegratedFissionSource = calculator::cell::IntegratedFissionSource<dim>;
  using ReturnType = eigenvalue::k_effective::UpdaterViaFissionSource;

  ReportBuildingComponant("K_Effective updater");
  std::unique_ptr<KEffectiveUpdater> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<ReturnType>(
          std::make_unique<AggregatedFissionSource>(
              std::make_unique<IntegratedFissionSource>(finite_element_ptr,
                                                        cross_sections_ptr),
              domain_ptr),
          2.0,
          10));

  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildMomentCalculator(MomentCalculatorImpl implementation)
-> std::unique_ptr<MomentCalculator> {
  return BuildMomentCalculator(nullptr, implementation);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildMomentCalculator(std::shared_ptr<QuadratureSet> quadrature_set_ptr,
                                                  FrameworkBuilder::MomentCalculatorImpl implementation)
-> std::unique_ptr<MomentCalculator> {
  ReportBuildingComponant("Moment calculator");
  std::unique_ptr<MomentCalculator> return_ptr = nullptr;

  try {
    return_ptr = std::move(quadrature::factory::MakeMomentCalculator<dim>(implementation, quadrature_set_ptr));

    if (implementation == MomentCalculatorImpl::kScalarMoment) {
      ReportBuildSuccess("(default) calculator for scalar solve");
    } else if (implementation == MomentCalculatorImpl::kZerothMomentOnly) {
      ReportBuildSuccess("(default) calculator for 0th moment only");
    } else {
      AssertThrow(false,
                  dealii::ExcMessage("Unsupported implementation of moment calculator specified in call to "
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
    double max_delta,
    int max_iterations) -> std::unique_ptr<MomentConvergenceChecker>{
  //TODO(Josh): Add option for using other than L1Norm
  ReportBuildingComponant("Moment convergence checker");

  using CheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  using FinalCheckerType = convergence::FinalCheckerOrN<system::moments::MomentVector,
                                                        convergence::moments::SingleMomentCheckerI>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(std::move(single_checker_ptr));
  return_ptr->SetMaxIterations(max_iterations);
  return return_ptr;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildMomentMapConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentMapConvergenceChecker> {
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
    std::unique_ptr<GroupSolveIteration> group_iteration_ptr,
    std::unique_ptr<ParameterConvergenceChecker> convergence_checker_ptr,
    const std::string&)
    -> std::unique_ptr<OuterIteration> {
  ReportBuildingComponant("Outer iteration");
  std::unique_ptr<OuterIteration> return_ptr = nullptr;
  using ReturnType = iteration::outer::OuterFixedSourceIteration;

  return_ptr = std::move(std::make_unique<ReturnType>(
      std::move(group_iteration_ptr),
      std::move(convergence_checker_ptr)));

  using ConvergenceDataPort = iteration::outer::data_names::ConvergenceStatusPort;
  using StatusPort =  iteration::outer::data_names::StatusPort;

  instrumentation::GetPort<ConvergenceDataPort>(*return_ptr)
      .AddInstrument(convergence_status_instrument_ptr_);

  instrumentation::GetPort<StatusPort>(*return_ptr)
      .AddInstrument(status_instrument_ptr_);

  ReportBuildSuccess(return_ptr->description());

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildOuterIteration(
    std::unique_ptr<GroupSolveIteration> group_solve_iteration_ptr,
    std::unique_ptr<ParameterConvergenceChecker> parameter_convergence_checker_ptr,
    std::unique_ptr<KEffectiveUpdater> k_effective_updater_ptr,
    const std::shared_ptr<FissionSourceUpdater>& fission_source_updater_ptr,
    const std::string& output_filename_base)
-> std::unique_ptr<OuterIteration> {
  std::unique_ptr<OuterIteration> return_ptr = nullptr;
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

  using InstrumentBuilder = instrumentation::builder::InstrumentBuilder;
  instrumentation::GetPort<ConvergenceDataPort>(*return_ptr)
      .AddInstrument(convergence_status_instrument_ptr_);
  instrumentation::GetPort<StatusPort>(*return_ptr)
      .AddInstrument(status_instrument_ptr_);
  instrumentation::GetPort<IterationErrorPort>(*return_ptr)
      .AddInstrument(Shared(InstrumentBuilder::BuildInstrument<std::pair<int,double>>(
          instrumentation::builder::InstrumentName::kIntDoublePairToFile,
          output_filename_base + "_iteration_error.csv")));

  validator_ptr_->AddPart(FrameworkPart::FissionSourceUpdate);
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildParameterConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<ParameterConvergenceChecker>{
  ReportBuildingComponant("Parameter (double) convergence checker");
  using CheckerType = convergence::parameters::SingleParameterChecker;
  using FinalCheckerType = convergence::FinalCheckerOrN<double, CheckerType>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));
  return_ptr->SetMaxIterations(max_iterations);

  return return_ptr;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildQuadratureSet(
    const problem::AngularQuadType quadrature_type,
    const FrameworkParameters::AngularQuadratureOrder order) -> std::shared_ptr<QuadratureSet> {
  ReportBuildingComponant("quadrature set");
  using QuadratureGenerator = quadrature::QuadratureGeneratorI<dim>;

  std::shared_ptr<QuadratureSet> return_ptr{ nullptr };
  std::shared_ptr<QuadratureGenerator> quadrature_generator_ptr{ nullptr };

  try {

    switch (quadrature_type) {
      case problem::AngularQuadType::kLevelSymmetricGaussian: {
        AssertThrow(dim == 3, dealii::ExcMessage("Error in BuildQuadratureSet LSGC only available for 3D"))
        quadrature_generator_ptr = quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
            order, quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
        break;
      }
      case problem::AngularQuadType::kGaussLegendre: {
        AssertThrow(dim == 1, dealii::ExcMessage("Error in BuildQuadratureSet GaussLegendre only available for 1D"))
        quadrature_generator_ptr =
            quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
                order, quadrature::AngularQuadratureSetType::kGaussLegendre);
        break;
      }
      default: {
        AssertThrow(false, dealii::ExcMessage("No supported quadratures for this dimension and transport model"))
        break;
      }
    }
    ReportBuildSuccess(quadrature_generator_ptr->description());

    return_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

    auto quadrature_points = quadrature::utility::GenerateAllPositiveX<dim>(quadrature_generator_ptr->GenerateSet());

    quadrature::factory::FillQuadratureSet<dim>(return_ptr.get(), quadrature_points);

    return return_ptr;
  } catch (...) {
    ReportBuildError();
    throw;
  }
}

template <int dim>
auto FrameworkBuilder<dim>::BuildSAAFFormulation(
    const std::shared_ptr<FiniteElement>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const std::shared_ptr<QuadratureSet>& quadrature_set_ptr,
    const formulation::SAAFFormulationImpl implementation) -> std::unique_ptr<SAAFFormulation> {
  ReportBuildingComponant("Building SAAF Formulation");
  std::unique_ptr<SAAFFormulation> return_ptr;

  if (implementation == formulation::SAAFFormulationImpl::kDefault) {
    using ReturnType = formulation::angular::SelfAdjointAngularFlux<dim>;
    return_ptr = std::move(std::make_unique<ReturnType>(finite_element_ptr, cross_sections_ptr, quadrature_set_ptr));
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSingleGroupSolver(const int max_iterations, const double convergence_tolerance)
-> std::unique_ptr<SingleGroupSolver> {
  using SolverName = solver::builder::SolverName;
  using SolverBuilder = solver::builder::SolverBuilder;

  ReportBuildingComponant("Single group solver");
  std::unique_ptr<SingleGroupSolver> return_ptr = nullptr;

  return_ptr = std::move(SolverBuilder::BuildSolver(SolverName::kDefaultGMRESGroupSolver, max_iterations,
                                                    convergence_tolerance));

  ReportBuildSuccess("Default implementation with GMRES");

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSystem(
    const int total_groups,
    const int total_angles,
    const Domain& domain,
    const std::size_t solution_size,
    bool is_eigenvalue_problem,
    bool need_rhs_boundary_condition) -> std::unique_ptr<System> {
  std::unique_ptr<System> return_ptr;

  ReportBuildingComponant("system");
  try {
    return_ptr = std::move(std::make_unique<System>());
    system_helper_.InitializeSystem(*return_ptr, total_groups, total_angles,
                             is_eigenvalue_problem, need_rhs_boundary_condition);
    system_helper_.SetUpSystemTerms(*return_ptr, domain);
    system_helper_.SetUpSystemMoments(*return_ptr, solution_size);
    ReportBuildSuccess("system");
  } catch (...) {
    ReportBuildError("system initialization error.");
    throw;
  }

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildStamper(const std::shared_ptr<Domain>& domain_ptr) -> std::unique_ptr<Stamper> {
  ReportBuildingComponant("Stamper");
  std::unique_ptr<Stamper> return_ptr = nullptr;

  return_ptr = std::move(std::make_unique<formulation::Stamper<dim>>(domain_ptr));
  ReportBuildSuccess(return_ptr->description());
  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildSubroutine(std::unique_ptr<FrameworkI> framework_ptr,
                                            const SubroutineName subroutine_name) -> std::unique_ptr<Subroutine> {
  std::unique_ptr<Subroutine> return_ptr{ nullptr };
  if (subroutine_name == SubroutineName::kGetScalarFluxFromFramework) {
    using ReturnType = iteration::subroutine::GetScalarFluxFromFramework;
    return_ptr = std::move(std::make_unique<ReturnType>(std::move(framework_ptr)));
  }
  return return_ptr;
}

template<int dim>
void FrameworkBuilder<dim>::Validate() const {
  validator_ptr_->ReportValidation();
}

template class FrameworkBuilder<1>;
template class FrameworkBuilder<2>;
template class FrameworkBuilder<3>;

} // namespace bart::framework::builder