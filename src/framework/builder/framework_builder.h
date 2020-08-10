#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_

#include <fstream>
#include <memory>
#include <data/cross_sections.h>
#include <deal.II/base/conditional_ostream.h>

#include "quadrature/quadrature_types.h"

#include "utility/has_description.h"

#include "framework/builder/framework_validator.h"
// Problem parameters
#include "problem/parameters_i.h"
#include "system/solution/solution_types.h"

// Interface classes built by this factory
#include "convergence/final_i.h"
#include "data/cross_sections.h"
#include "domain/definition_i.h"
#include "domain/finite_element/finite_element_i.h"
#include "eigenvalue/k_effective/k_effective_updater_i.h"
#include "formulation/stamper_i.h"
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.h"
#include "formulation/updater/fission_source_updater_i.h"
#include "formulation/updater/fixed_updater_i.h"
#include "formulation/updater/scattering_source_updater_i.h"
#include "formulation/updater/boundary_conditions_updater_i.h"
#include "framework/framework_i.h"
#include "instrumentation/instrument_i.h"
#include "iteration/group/group_solve_iteration_i.h"
#include "iteration/initializer/initializer_i.h"
#include "iteration/outer/outer_iteration_i.h"
#include "instrumentation/port.h"
#include "quadrature/quadrature_set_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "solver/group/single_group_solver_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/system.h"

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

  using AngularFluxStorage = system::solution::EnergyGroupToAngularSolutionPtrMap;

  using BoundaryConditionsUpdaterType = formulation::updater::BoundaryConditionsUpdaterI;
  using ConvergenceInstrumentType = instrumentation::InstrumentI<convergence::Status>;
  using CrossSectionType = data::CrossSections;
  using DiffusionFormulationType = formulation::scalar::DiffusionI<dim>;
  using DomainType = domain::DefinitionI<dim>;
  using FiniteElementType = domain::finite_element::FiniteElementI<dim>;
  using FissionSourceUpdaterType = formulation::updater::FissionSourceUpdaterI;
  using FixedUpdaterType = formulation::updater::FixedUpdaterI;
  using FrameworkType = framework::FrameworkI;
  using GroupSolutionType = system::solution::MPIGroupAngularSolutionI;
  using GroupSolveIterationType = iteration::group::GroupSolveIterationI;
  using InitializerType = iteration::initializer::InitializerI;
  using KEffectiveUpdaterType = eigenvalue::k_effective::K_EffectiveUpdaterI;
  using MomentCalculatorType = quadrature::calculators::SphericalHarmonicMomentsI;
  using MomentConvergenceCheckerType = convergence::FinalI<system::moments::MomentVector>;
  using MomentMapConvergenceCheckerType = convergence::FinalI<const system::moments::MomentsMap>;
  using OuterIterationType = iteration::outer::OuterIterationI;
  using ParameterConvergenceCheckerType = convergence::FinalI<double>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using ScatteringSourceUpdaterType = formulation::updater::ScatteringSourceUpdaterI;
  using SingleGroupSolverType = solver::group::SingleGroupSolverI;
  using StamperType = formulation::StamperI<dim>;
  using SystemType = system::System;

  struct UpdaterPointers {
    std::shared_ptr<BoundaryConditionsUpdaterType> boundary_conditions_updater_ptr = nullptr;
    std::shared_ptr<FissionSourceUpdaterType> fission_source_updater_ptr = nullptr;
    std::shared_ptr<FixedUpdaterType> fixed_updater_ptr = nullptr;
    std::shared_ptr<ScatteringSourceUpdaterType> scattering_source_updater_ptr = nullptr;
  };

  FrameworkBuilder(std::shared_ptr<FrameworkReporterType> reporter_ptr)
  : reporter_ptr_(reporter_ptr) {}
  ~FrameworkBuilder() = default;

  std::unique_ptr<FrameworkType> BuildFramework(std::string name, ParametersType&);

  std::unique_ptr<ConvergenceInstrumentType> BuildConvergenceInstrument();
  std::unique_ptr<CrossSectionType> BuildCrossSections(ParametersType);
  std::unique_ptr<DiffusionFormulationType> BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElementType>&,
      const std::shared_ptr<data::CrossSections>&,
      const formulation::DiffusionFormulationImpl implementation = formulation::DiffusionFormulationImpl::kDefault);
  std::unique_ptr<DomainType> BuildDomain(
      ParametersType, const std::shared_ptr<FiniteElementType>&,
      std::string material_mapping);
  std::unique_ptr<FiniteElementType> BuildFiniteElement(ParametersType);
  UpdaterPointers BuildUpdaterPointers(
      std::unique_ptr<DiffusionFormulationType>,
      std::unique_ptr<StamperType>,
      const std::map<problem::Boundary, bool>& reflective_boundaries);
  UpdaterPointers BuildUpdaterPointers(
      std::unique_ptr<SAAFFormulationType>,
      std::unique_ptr<StamperType>,
      const std::shared_ptr<QuadratureSetType>&);
  UpdaterPointers BuildUpdaterPointers(
      std::unique_ptr<SAAFFormulationType>,
      std::unique_ptr<StamperType>,
      const std::shared_ptr<QuadratureSetType>&,
      const std::map<problem::Boundary, bool>& reflective_boundaries,
      const AngularFluxStorage&);
  std::unique_ptr<GroupSolveIterationType> BuildGroupSolveIteration(
      std::unique_ptr<SingleGroupSolverType>,
      std::unique_ptr<MomentConvergenceCheckerType>,
      std::unique_ptr<MomentCalculatorType>,
      const std::shared_ptr<GroupSolutionType>&,
      const UpdaterPointers& updater_ptrs,
      std::unique_ptr<MomentMapConvergenceCheckerType> moment_map_convergence_checker_ptr);
  std::unique_ptr<InitializerType> BuildInitializer(
      const std::shared_ptr<formulation::updater::FixedUpdaterI>&,
      const int total_groups, const int total_angles);
  std::unique_ptr<GroupSolutionType> BuildGroupSolution(const int n_angles);
  std::unique_ptr<KEffectiveUpdaterType> BuildKEffectiveUpdater(
      const std::shared_ptr<FiniteElementType>&,
      const std::shared_ptr<CrossSectionType>&,
      const std::shared_ptr<DomainType>&);
  std::unique_ptr<MomentCalculatorType> BuildMomentCalculator(
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kScalarMoment);
  std::unique_ptr<MomentCalculatorType> BuildMomentCalculator(
      std::shared_ptr<QuadratureSetType>,
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kZerothMomentOnly);
  std::unique_ptr<MomentConvergenceCheckerType> BuildMomentConvergenceChecker(
      double max_delta, int max_iterations);
  std::unique_ptr<MomentMapConvergenceCheckerType> BuildMomentMapConvergenceChecker(
      double max_delta, int max_iterations);
  std::unique_ptr<OuterIterationType> BuildOuterIteration(
      std::unique_ptr<GroupSolveIterationType>,
      std::unique_ptr<ParameterConvergenceCheckerType>);
  std::unique_ptr<OuterIterationType> BuildOuterIteration(
      std::unique_ptr<GroupSolveIterationType>,
      std::unique_ptr<ParameterConvergenceCheckerType>,
      std::unique_ptr<KEffectiveUpdaterType>,
      const std::shared_ptr<FissionSourceUpdaterType>&);
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
  std::unique_ptr<SystemType> BuildSystem(const int n_groups, const int n_angles,
                                          const DomainType& domain,
                                          const std::size_t solution_size,
                                          bool is_eigenvalue_problem = true,
                                          bool need_rhs_boundary_condition = false);

  FrameworkReporterType* reporter_ptr() { return reporter_ptr_.get(); }

 private:
  void ReportBuildingComponant(std::string componant) {
    if (!build_report_closed_) {
      *reporter_ptr_ << "\n" << Color::Reset;
    }
    *reporter_ptr_ << "\tBuilding " << componant << ": ";
    build_report_closed_ = false;
  }

  void Report(const std::string to_report, const Color color = Color::Reset) const {
    *reporter_ptr_ << color << to_report << Color::Reset;
  }

  void ReportBuildSuccess(std::string description = "") {
    *reporter_ptr_ << Color::Green << "Built " << description << Color::Reset
                   << "\n";
    build_report_closed_ = true;
  }

  void ReportBuildError(std::string description = "") {
    *reporter_ptr_ << Color::Red << "Error: " << description << Color::Reset
                   << "\n";
    build_report_closed_ = true;
  }

  void Validate() const;

  std::shared_ptr<FrameworkReporterType> reporter_ptr_;

  template <typename T>
  inline std::shared_ptr<T> Shared(std::unique_ptr<T> to_convert_ptr) {
    return to_convert_ptr;
  }

  std::string ReadMappingFile(std::string filename);

  FrameworkValidator validator_;
  bool build_report_closed_ = true;
};

} // namespace builder

} // namespace framework

} // namespace bart

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_H_
