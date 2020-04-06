#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include <fstream>
#include <memory>
#include <data/cross_sections.h>
#include <deal.II/base/conditional_ostream.h>

#include "quadrature/quadrature_types.h"

// Problem parameters
#include "problem/parameters_i.h"

// Interface classes built by this factory
#include "convergence/reporter/mpi_i.h"
#include "convergence/final_i.h"
#include "data/cross_sections.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/stamper_i.h"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "iteration/initializer/initializer_i.h"
#include "quadrature/quadrature_set_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "solver/group/single_group_solver_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"

// Dependency clases
#include "formulation/updater/fixed_updater_i.h"
#include "utility/reporter/basic_reporter_i.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
class FrameworkBuilder {
 public:
  using FrameworkReporterType = utility::reporter::BasicReporterI;
  using ParametersType = const problem::ParametersI&;
  using Color = utility::reporter::Color;
  using MomentCalculatorImpl = quadrature::MomentCalculatorImpl;

  using CrossSectionType = data::CrossSections;
  using DiffusionFormulationType = formulation::scalar::DiffusionI<dim>;
  using DomainType = domain::DefinitionI<dim>;
  using FiniteElementType = domain::finite_element::FiniteElementI<dim>;
  using FixedUpdaterType = formulation::updater::FixedUpdaterI;
  using GroupSolutionType = system::solution::MPIGroupAngularSolutionI;
  using InitializerType = iteration::initializer::InitializerI;
  using MomentCalculatorType = quadrature::calculators::SphericalHarmonicMomentsI;
  using MomentConvergenceCheckerType = convergence::FinalI<system::moments::MomentVector>;
  using ParameterConvergenceCheckerType = convergence::FinalI<double>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using ReporterType = convergence::reporter::MpiI;
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using SingleGroupSolverType = solver::group::SingleGroupSolverI;
  using StamperType = formulation::StamperI<dim>;


  FrameworkBuilder(std::shared_ptr<FrameworkReporterType> reporter_ptr)
  : reporter_ptr_(reporter_ptr) {}
  ~FrameworkBuilder() = default;

  void BuildFramework(std::string name, ParametersType&);

  std::unique_ptr<ReporterType> BuildConvergenceReporter();
  std::unique_ptr<CrossSectionType> BuildCrossSections(ParametersType);
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
  std::unique_ptr<GroupSolutionType> BuildGroupSolution(const int n_angles);
  std::unique_ptr<MomentCalculatorType> BuildMomentCalculator(
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kScalarMoment);
  std::unique_ptr<MomentCalculatorType> BuildMomentCalculator(
      std::shared_ptr<QuadratureSetType>,
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kZerothMomentOnly);
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

  FrameworkReporterType* reporter_ptr() { return reporter_ptr_.get(); }

 private:
  void ReportBuildingComponant(std::string componant) {
    reporter_ptr_->Report("\tBuilding " + componant + ": "); }
  void ReportBuilt(std::string description) {
    reporter_ptr_->Report("Built " + description + "\n", Color::Green); }
  void ReportError() {
    reporter_ptr_->Report("Error\n", Color::Red); }

  std::shared_ptr<FrameworkReporterType> reporter_ptr_;
  template <typename T>
  inline std::shared_ptr<T> Shared(std::unique_ptr<T> to_convert_ptr) {
    return std::move(to_convert_ptr);
  }
  std::string ReadMappingFile(std::string filename);
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
