#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_

#include <memory>
#include <string>

#include "convergence/final_i.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "framework/framework_i.hpp"
#include "framework/framework_parameters.hpp"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.h"
#include "formulation/updater/boundary_conditions_updater_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/stamper_i.h"
#include "iteration/initializer/initializer_i.h"
#include "iteration/group/group_solve_iteration_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "quadrature/quadrature_set_i.h"
#include "problem/parameter_types.h"
#include "solver/group/single_group_solver_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/solution/solution_types.h"

namespace bart::framework::builder {

template <int dim>
class FrameworkBuilderI {
 public:
  // Classes built by member functions
  using DiffusionFormulation = typename formulation::scalar::DiffusionI<dim>;
  using Domain = typename domain::DefinitionI<dim>;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;
  using FrameworkI = framework::FrameworkI;
  using GroupSolution = system::solution::MPIGroupAngularSolutionI;
  using GroupSolveIteration = iteration::group::GroupSolveIterationI;
  using Initializer = iteration::initializer::InitializerI;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsI;
  using MomentConvergenceChecker = convergence::FinalI<system::moments::MomentVector>;
  using MomentMapConvergenceChecker = convergence::FinalI<const system::moments::MomentsMap>;
  using QuadratureSet = typename quadrature::QuadratureSetI<dim>;
  using SAAFFormulation = typename formulation::angular::SelfAdjointAngularFluxI<dim>;
  using SingleGroupSolver = solver::group::SingleGroupSolverI;
  using Stamper = formulation::StamperI<dim>;

  // Other types
  using AngularFluxStorage = system::solution::EnergyGroupToAngularSolutionPtrMap;

  // Implementation specifiers
  using DiffusionFormulationImpl = formulation::DiffusionFormulationImpl;
  using MomentCalculatorImpl = quadrature::MomentCalculatorImpl;

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

  virtual auto BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::CrossSections>&,
      const DiffusionFormulationImpl) -> std::unique_ptr<DiffusionFormulation> = 0;
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
  virtual auto BuildMomentCalculator(MomentCalculatorImpl) -> std::unique_ptr<MomentCalculator> = 0;
  virtual auto BuildMomentCalculator(std::shared_ptr<QuadratureSet>,
                                     MomentCalculatorImpl) -> std::unique_ptr<MomentCalculator> = 0;
  virtual auto BuildMomentConvergenceChecker(double max_delta,
                                             int max_iterations) -> std::unique_ptr<MomentConvergenceChecker> = 0;
  virtual auto BuildMomentMapConvergenceChecker(double max_delta,
                                                int max_iterations) -> std::unique_ptr<MomentMapConvergenceChecker> = 0;
  virtual auto BuildQuadratureSet(
      const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder) -> std::shared_ptr<QuadratureSet> = 0;
  virtual auto BuildSAAFFormulation(const std::shared_ptr<FiniteElement>&,
                                    const std::shared_ptr<data::CrossSections>&,
                                    const std::shared_ptr<QuadratureSet>&,
                                    const formulation::SAAFFormulationImpl) -> std::unique_ptr<SAAFFormulation> = 0;
  virtual auto BuildSingleGroupSolver(const int max_iterations,
                                      const double convergence_tolerance) -> std::unique_ptr<SingleGroupSolver> = 0;
  virtual auto BuildStamper(const std::shared_ptr<Domain>&) -> std::unique_ptr<Stamper> = 0;
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
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_I_HPP_
