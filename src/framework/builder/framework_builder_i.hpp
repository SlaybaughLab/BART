#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_

#include <memory>
#include <string>

#include "data/cross_sections/material_cross_sections.hpp"
#include "convergence/status.hpp"
#include "convergence/iteration_completion_checker_i.hpp"
#include "domain/domain_i.hpp"
#include "domain/finite_element/finite_element_i.hpp"
#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"
#include "framework/framework_i.hpp"
#include "framework/framework_parameters.hpp"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.hpp"
#include "formulation/scalar/drift_diffusion_i.hpp"
#include "formulation/updater/boundary_conditions_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/stamper_i.hpp"
#include "framework/builder/framework_validator.hpp"
#include "instrumentation/instrument_i.h"
#include "iteration/initializer/initializer_i.h"
#include "iteration/group/group_solve_iteration_i.hpp"
#include "iteration/outer/outer_iteration_i.hpp"
#include "iteration/subroutine/subroutine_i.hpp"
#include "iteration/initializer/factory.hpp"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "quadrature/calculators/angular_flux_integrator_i.hpp"
#include "quadrature/quadrature_set_i.hpp"
#include "problem/parameter_types.hpp"
#include "solver/group/single_group_solver_i.h"
#include "system/moments/spherical_harmonic_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/solution/solution_types.h"
#include "system/system.hpp"
#include "utility/colors.hpp"

namespace bart::framework::builder {

template <int dim>
class FrameworkBuilderI {
 public:
  // Classes built by member functions
  using AngularFluxIntegrator = quadrature::calculators::AngularFluxIntegratorI;
  using CrossSections = data::cross_sections::MaterialCrossSections;
  using DiffusionFormulation = typename formulation::scalar::DiffusionI<dim>;
  using DriftDiffusionFormulation = typename formulation::scalar::DriftDiffusionI<dim>;
  using Domain = typename domain::DomainI<dim>;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using FrameworkI = framework::FrameworkI;
  using GroupSolution = system::solution::MPIGroupAngularSolutionI;
  using GroupSolveIteration = iteration::group::GroupSolveIterationI;
  using Initializer = iteration::initializer::InitializerI;
  using KEffectiveUpdater = eigenvalue::k_eigenvalue::K_EigenvalueCalculatorI;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsI;
  using MomentConvergenceChecker = convergence::IterationCompletionCheckerI<system::moments::MomentVector>;
  using MomentMapConvergenceChecker = convergence::IterationCompletionCheckerI<system::moments::MomentsMap>;
  using OuterIteration = iteration::outer::OuterIterationI;
  using ParameterConvergenceChecker = convergence::IterationCompletionCheckerI<double>;
  using QuadratureSet = typename quadrature::QuadratureSetI<dim>;
  using SAAFFormulation = typename formulation::angular::SelfAdjointAngularFluxI<dim>;
  using SphericalHarmonicMoments = system::moments::SphericalHarmonicI;
  using SingleGroupSolver = solver::group::SingleGroupSolverI;
  using Stamper = formulation::StamperI<dim>;
  using Subroutine = iteration::subroutine::SubroutineI;
  using System = system::System;
  using Validator = framework::builder::FrameworkValidatorI;

  // Instrumentation
  using ColorStatusPair = std::pair<std::string, utility::Color>;
  using ColorStatusInstrument = instrumentation::InstrumentI<ColorStatusPair>;
  using ConvergenceInstrument = instrumentation::InstrumentI<convergence::Status>;
  using StatusInstrument = instrumentation::InstrumentI<std::string>;

  // Other types
  using AngularFluxStorage = system::solution::EnergyGroupToAngularSolutionPtrMap;

  // Implementation specifiers
  using DiffusionFormulationImpl = formulation::DiffusionFormulationImpl;
  using InitializerName = iteration::initializer::InitializerName;
  using MomentCalculatorImpl = quadrature::MomentCalculatorImpl;
  using SubroutineName = iteration::subroutine::SubroutineName;

  // Updater Pointers
  using BoundaryConditionsUpdater = formulation::updater::BoundaryConditionsUpdaterI;
  using FissionSourceUpdater = formulation::updater::FissionSourceUpdaterI;
  using FixedTermUpdater = formulation::updater::FixedUpdaterI;
  using ScatteringSourceUpdater = formulation::updater::ScatteringSourceUpdaterI;

  struct UpdaterPointers {
    std::shared_ptr<BoundaryConditionsUpdater> boundary_conditions_updater_ptr{ nullptr };
    std::shared_ptr<FissionSourceUpdater> fission_source_updater_ptr{ nullptr };
    std::shared_ptr<FixedTermUpdater> fixed_updater_ptr{ nullptr };
    std::shared_ptr<ScatteringSourceUpdater> scattering_source_updater_ptr{ nullptr };
  };

  virtual ~FrameworkBuilderI() = default;

  virtual auto BuildAngularFluxIntegrator(
      const std::shared_ptr<QuadratureSet>) -> std::unique_ptr<AngularFluxIntegrator> = 0;
  virtual auto BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
      const DiffusionFormulationImpl) -> std::unique_ptr<DiffusionFormulation> = 0;
  virtual auto BuildDriftDiffusionFormulation(
      const std::shared_ptr<AngularFluxIntegrator>&,
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::cross_sections::MaterialCrossSections>&) -> std::unique_ptr<DriftDiffusionFormulation> = 0;
  virtual auto BuildDomain(const FrameworkParameters::DomainSize,
                           const FrameworkParameters::NumberOfCells,
                           const std::shared_ptr<FiniteElement>&,
                           const std::string material_mapping) -> std::unique_ptr<Domain> = 0;
  virtual auto BuildFiniteElement(
      const problem::CellFiniteElementType,
      const problem::DiscretizationType,
      const FrameworkParameters::PolynomialDegree) -> std::unique_ptr<FiniteElement> = 0;
  virtual auto BuildGroupSolution(const int n_angles) -> std::unique_ptr<GroupSolution> = 0;
  virtual auto BuildGroupSolveIteration(
      std::unique_ptr<SingleGroupSolver>,
      std::unique_ptr<MomentConvergenceChecker>,
      std::unique_ptr<MomentCalculator>,
      const std::shared_ptr<GroupSolution>&,
      const UpdaterPointers& updater_ptrs,
      std::unique_ptr<MomentMapConvergenceChecker>) -> std::unique_ptr<GroupSolveIteration> = 0;
  virtual auto BuildInitializer(const std::shared_ptr<FixedTermUpdater>&,
                                const int total_groups,
                                const int total_angles) -> std::unique_ptr<Initializer> = 0;
  virtual auto BuildInitializer(const std::shared_ptr<FixedTermUpdater>&,
                                const int total_groups,
                                const int total_angles,
                                const InitializerName) -> std::unique_ptr<Initializer> = 0;
  virtual auto BuildKEffectiveUpdater() -> std::unique_ptr<KEffectiveUpdater> = 0;
  virtual auto BuildKEffectiveUpdater(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<CrossSections>&,
      const std::shared_ptr<Domain>&) -> std::unique_ptr<KEffectiveUpdater> = 0;
  virtual auto BuildMomentCalculator(MomentCalculatorImpl) -> std::unique_ptr<MomentCalculator> = 0;
  virtual auto BuildMomentCalculator(std::shared_ptr<QuadratureSet>,
                                     MomentCalculatorImpl) -> std::unique_ptr<MomentCalculator> = 0;
  virtual auto BuildMomentConvergenceChecker(double max_delta,
                                             int max_iterations) -> std::unique_ptr<MomentConvergenceChecker> = 0;
  virtual auto BuildMomentMapConvergenceChecker(double max_delta,
                                                int max_iterations) -> std::unique_ptr<MomentMapConvergenceChecker> = 0;
  virtual auto BuildOuterIteration(std::unique_ptr<GroupSolveIteration>,
                                   std::unique_ptr<ParameterConvergenceChecker>,
                                   const std::string& output_filename_base) -> std::unique_ptr<OuterIteration> = 0;
  virtual auto BuildOuterIteration(
      std::unique_ptr<GroupSolveIteration>,
      std::unique_ptr<ParameterConvergenceChecker>,
      std::unique_ptr<KEffectiveUpdater>,
      const std::shared_ptr<FissionSourceUpdater>&,
      const std::string& output_filename_base) -> std::unique_ptr<OuterIteration> = 0;
  virtual auto BuildParameterConvergenceChecker(double max_delta,
                                                int max_iterations) -> std::unique_ptr<ParameterConvergenceChecker> = 0;
  virtual auto BuildQuadratureSet(
      const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder) -> std::shared_ptr<QuadratureSet> = 0;
  virtual auto BuildSAAFFormulation(const std::shared_ptr<FiniteElement>&,
                                    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
                                    const std::shared_ptr<QuadratureSet>&,
                                    const formulation::SAAFFormulationImpl) -> std::unique_ptr<SAAFFormulation> = 0;
  virtual auto BuildSingleGroupSolver(const int max_iterations,
                                      const double convergence_tolerance) -> std::unique_ptr<SingleGroupSolver> = 0;
  virtual auto BuildStamper(const std::shared_ptr<Domain>&) -> std::unique_ptr<Stamper> = 0;
  virtual auto BuildSubroutine(std::unique_ptr<FrameworkI>, const SubroutineName) -> std::unique_ptr<Subroutine> = 0;
  virtual auto BuildSystem(const int n_groups,
                           const int n_angles,
                           const Domain& domain,
                           const std::size_t solution_size,
                           bool is_eigenvalue_problem,
                           bool need_rhs_boundary_condition) -> std::unique_ptr<System> = 0;

  virtual auto BuildUpdaterPointers(std::unique_ptr<DiffusionFormulation>,
                                    std::unique_ptr<DriftDiffusionFormulation>,
                                    std::shared_ptr<Stamper>,
                                    std::shared_ptr<AngularFluxIntegrator>,
                                    std::shared_ptr<SphericalHarmonicMoments>,
                                    AngularFluxStorage&,
                                    const std::map<problem::Boundary, bool>&) -> UpdaterPointers = 0;
  virtual auto BuildUpdaterPointers(std::unique_ptr<DiffusionFormulation>,
                                    std::unique_ptr<Stamper>,
                                    const std::map<problem::Boundary, bool>&) -> UpdaterPointers = 0;
  virtual auto BuildUpdaterPointers(std::unique_ptr<SAAFFormulation>,
                                    std::unique_ptr<Stamper>,
                                    const std::shared_ptr<QuadratureSet>&) -> UpdaterPointers = 0;
  virtual auto BuildUpdaterPointers(std::unique_ptr<SAAFFormulation>,
                                    std::unique_ptr<Stamper>,
                                    const std::shared_ptr<QuadratureSet>&,
                                    const std::map<problem::Boundary, bool>& reflective_boundaries,
                                    const AngularFluxStorage&) -> UpdaterPointers = 0;

  // Instrumentation getters and setters
  virtual auto set_color_status_instrument_ptr(
      const std::shared_ptr<ColorStatusInstrument>&) -> FrameworkBuilderI<dim>& = 0;
  virtual auto set_convergence_status_instrument_ptr(
      const std::shared_ptr<ConvergenceInstrument>&) -> FrameworkBuilderI<dim>& = 0;
  virtual auto set_status_instrument_ptr(
      const std::shared_ptr<StatusInstrument>&) -> FrameworkBuilderI<dim>& = 0;
  virtual auto color_status_instrument_ptr() const -> std::shared_ptr<ColorStatusInstrument> = 0;
  virtual auto convergence_status_instrument_ptr() const -> std::shared_ptr<ConvergenceInstrument> = 0;
  virtual auto status_instrument_ptr() const -> std::shared_ptr<StatusInstrument> = 0;

  // Access the internal validator object
  virtual auto validator_ptr() -> Validator* = 0;
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
