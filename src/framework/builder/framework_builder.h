#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include <memory>
#include <data/cross_sections.h>

// Problem parameters
#include "problem/parameters_i.h"

// Interface classes built by this factory
#include "convergence/reporter/mpi_i.h"
#include "convergence/final_i.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/stamper_i.h"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "iteration/initializer/initializer_i.h"
#include "quadrature/quadrature_set_i.h"
#include "solver/group/single_group_solver_i.h"

// Dependency clases
#include "formulation/updater/fixed_updater_i.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class FrameworkBuilder {
 public:
  FrameworkBuilder() = default;
  ~FrameworkBuilder() = default;

  using ParametersType = const problem::ParametersI&;

  using DiffusionFormulationType = formulation::scalar::DiffusionI<dim>;
  using DomainType = domain::DefinitionI<dim>;
  using FiniteElementType = domain::finite_element::FiniteElementI<dim>;
  using FixedUpdaterType = formulation::updater::FixedUpdaterI;
  using InitializerType = iteration::initializer::InitializerI;
  using MomentConvergenceCheckerType = convergence::FinalI<system::moments::MomentVector>;
  using ParameterConvergenceCheckerType = convergence::FinalI<double>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using ReporterType = convergence::reporter::MpiI;
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using SingleGroupSolverType = solver::group::SingleGroupSolverI;
  using StamperType = formulation::StamperI<dim>;

  std::unique_ptr<ReporterType> BuildConvergenceReporter();
  std::unique_ptr<DiffusionFormulationType> BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElementType>&,
      const std::shared_ptr<data::CrossSections>&,
      const formulation::DiffusionFormulationImpl implementation = formulation::DiffusionFormulationImpl::kDefault);
  std::unique_ptr<DomainType> BuildDomain(
      ParametersType, const std::shared_ptr<FiniteElementType>&,
      std::string material_mapping);
  std::unique_ptr<FiniteElementType> BuildFiniteElement(ParametersType);
  std::unique_ptr<FixedUpdaterType> BuildFixedUpdater(
      std::unique_ptr<DiffusionFormulationType>,
      std::unique_ptr<StamperType>);
  std::unique_ptr<FixedUpdaterType> BuildFixedUpdater(
      std::unique_ptr<SAAFFormulationType>,
      std::unique_ptr<StamperType>,
      const std::shared_ptr<QuadratureSetType>&);
  std::unique_ptr<InitializerType> BuildInitializer(
      const std::shared_ptr<formulation::updater::FixedUpdaterI>&,
      const int total_groups, const int total_angles);
  std::unique_ptr<MomentConvergenceCheckerType> BuildMomentConvergenceChecker(
      double max_delta, int max_iterations);
  std::unique_ptr<ParameterConvergenceCheckerType> BuildParameterConvergenceChecker(
      double max_delta, int max_iterations);
  std::shared_ptr<QuadratureSetType> BuildQuadratureSet(ParametersType);
  std::unique_ptr<SAAFFormulationType> BuildSAAFFormulation(
      const std::shared_ptr<FiniteElementType>&,
      const std::shared_ptr<data::CrossSections>&,
      const std::shared_ptr<QuadratureSetType>&,
      const formulation::SAAFFormulationImpl implementation = formulation::SAAFFormulationImpl::kDefault);
  std::unique_ptr<SingleGroupSolverType> BuildSingleGroupSolver(
      const int max_iterations = 1000,
      const double convergence_tolerance = 1e-10);
  std::unique_ptr<StamperType> BuildStamper(const std::shared_ptr<DomainType>&);

};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
