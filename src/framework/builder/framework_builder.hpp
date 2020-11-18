#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_

#include "framework/builder/framework_builder_i.hpp"

#include <fstream>
#include <memory>
#include <data/cross_sections.h>
#include <deal.II/base/conditional_ostream.h>

#include "quadrature/quadrature_types.h"

#include "utility/has_description.h"
#include "system/system_helper.hpp"

#include "framework/framework_parameters.hpp"
#include "framework/builder/framework_validator.hpp"
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
#include "framework/framework_i.hpp"
#include "iteration/group/group_solve_iteration_i.h"
#include "iteration/outer/outer_iteration_i.hpp"
#include "instrumentation/port.hpp"
#include "instrumentation/instrument_i.h"
#include "quadrature/quadrature_set_i.h"
#include "solver/group/single_group_solver_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/system.h"
#include "system/moments/spherical_harmonic_i.h"

// Dependency clases
#include "formulation/updater/fixed_updater_i.h"
#include "utility/colors.hpp"
#include "instrumentation/port.hpp"

namespace bart::framework::builder {

namespace data_port {
struct BuilderStatus;
using StatusDataPort = instrumentation::Port<std::pair<std::string, utility::Color>, BuilderStatus>;
}

template <int dim>
class FrameworkBuilder : public data_port::StatusDataPort, public FrameworkBuilderI<dim> {
 public:
  // New using types from refactor
  using typename FrameworkBuilderI<dim>::DiffusionFormulation;
  using typename FrameworkBuilderI<dim>::Domain;
  using typename FrameworkBuilderI<dim>::FiniteElement;
  using typename FrameworkBuilderI<dim>::FrameworkI;
  using typename FrameworkBuilderI<dim>::Initializer;
  using typename FrameworkBuilderI<dim>::MomentCalculator;
  using typename FrameworkBuilderI<dim>::QuadratureSet;
  using typename FrameworkBuilderI<dim>::Stamper;
  using typename FrameworkBuilderI<dim>::SAAFFormulation;

  using typename FrameworkBuilderI<dim>::UpdaterPointers;
  using typename FrameworkBuilderI<dim>::BoundaryConditionsUpdater;
  using typename FrameworkBuilderI<dim>::FissionSourceUpdater;
  using typename FrameworkBuilderI<dim>::FixedTermUpdater;
  using typename FrameworkBuilderI<dim>::ScatteringSourceUpdater;

  using typename FrameworkBuilderI<dim>::DiffusionFormulationImpl;
  using typename FrameworkBuilderI<dim>::MomentCalculatorImpl;

  // TODO: Remove old types as they are unneeded
  using ParametersType = const problem::ParametersI&;
  using Color = utility::Color;

  using AngularFluxStorage = system::solution::EnergyGroupToAngularSolutionPtrMap;


  using CrossSectionType = data::CrossSections;
  using DomainType = domain::DefinitionI<dim>;
  using FiniteElementType = domain::finite_element::FiniteElementI<dim>;
  using FrameworkType = framework::FrameworkI;
  using GroupSolutionType = system::solution::MPIGroupAngularSolutionI;
  using GroupSolveIterationType = iteration::group::GroupSolveIterationI;
  using KEffectiveUpdaterType = eigenvalue::k_effective::K_EffectiveUpdaterI;
  using MomentConvergenceCheckerType = convergence::FinalI<system::moments::MomentVector>;
  using MomentMapConvergenceCheckerType = convergence::FinalI<const system::moments::MomentsMap>;
  using OuterIterationType = iteration::outer::OuterIterationI;
  using ParameterConvergenceCheckerType = convergence::FinalI<double>;
  using QuadratureSetType = quadrature::QuadratureSetI<dim>;
  using SAAFFormulationType = formulation::angular::SelfAdjointAngularFluxI<dim>;
  using SingleGroupSolverType = solver::group::SingleGroupSolverI;
  using StamperType = formulation::StamperI<dim>;
  using SystemType = system::System;

  using ColorStatusPair = std::pair<std::string, utility::Color>;
  // Instrumentation
  using ColorStatusInstrument = instrumentation::InstrumentI<ColorStatusPair>;
  using ConvergenceInstrument = instrumentation::InstrumentI<convergence::Status>;
  using StatusInstrument = instrumentation::InstrumentI<std::string>;

  FrameworkBuilder() = default;
  ~FrameworkBuilder() = default;

  std::unique_ptr<FrameworkType> BuildFramework(std::string name, ParametersType&);
  std::unique_ptr<FrameworkType> BuildFramework(std::string name, ParametersType&,
                                                system::moments::SphericalHarmonicI*);

  [[nodiscard]] auto BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::CrossSections>&,
      const DiffusionFormulationImpl implementation = DiffusionFormulationImpl::kDefault)
  -> std::unique_ptr<DiffusionFormulation> override;
  [[nodiscard]] auto BuildDomain(const FrameworkParameters::DomainSize,
                                 const FrameworkParameters::NumberOfCells,
                                 const std::shared_ptr<FiniteElement>&,
                                 const std::string material_mapping) -> std::unique_ptr<Domain> override;
  [[nodiscard]] auto BuildFiniteElement(
      const problem::CellFiniteElementType finite_element_type,
      const problem::DiscretizationType discretization_type,
      const FrameworkParameters::PolynomialDegree polynomial_degree) -> std::unique_ptr<FiniteElement> override;
  [[nodiscard]] auto BuildMomentCalculator(
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kScalarMoment)
  -> std::unique_ptr<MomentCalculator> override;
  [[nodiscard]] auto BuildMomentCalculator(
      std::shared_ptr<QuadratureSet>,
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kZerothMomentOnly)
  -> std::unique_ptr<MomentCalculator> override;
  [[nodiscard]] auto BuildQuadratureSet(
      const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder) -> std::shared_ptr<QuadratureSet> override;
  [[nodiscard]] auto BuildSAAFFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::CrossSections>&,
      const std::shared_ptr<QuadratureSet>&,
      const formulation::SAAFFormulationImpl implementation = formulation::SAAFFormulationImpl::kDefault)
  -> std::unique_ptr<SAAFFormulation> override;
  [[nodiscard]] auto BuildStamper(const std::shared_ptr<Domain>&) -> std::unique_ptr<Stamper> override;
  [[nodiscard]] auto BuildUpdaterPointers(
      std::unique_ptr<DiffusionFormulation>,
      std::unique_ptr<Stamper>,
      const std::map<problem::Boundary, bool>& reflective_boundaries) -> UpdaterPointers override;
  [[nodiscard]] auto BuildUpdaterPointers(std::unique_ptr<SAAFFormulation>,
                                          std::unique_ptr<Stamper>,
                                          const std::shared_ptr<QuadratureSet>&) -> UpdaterPointers override;
  [[nodiscard]] auto BuildUpdaterPointers(std::unique_ptr<SAAFFormulation>,
                                          std::unique_ptr<Stamper>,
                                          const std::shared_ptr<QuadratureSet>&,
                                          const std::map<problem::Boundary, bool>& reflective_boundaries,
                                          const AngularFluxStorage&) -> UpdaterPointers override;

  std::unique_ptr<CrossSectionType> BuildCrossSections(ParametersType);

  std::unique_ptr<DomainType> BuildDomain(
      ParametersType, const std::shared_ptr<FiniteElementType>&,
      std::string material_mapping);
  std::unique_ptr<FiniteElementType> BuildFiniteElement(ParametersType);
  std::unique_ptr<GroupSolveIterationType> BuildGroupSolveIteration(
      std::unique_ptr<SingleGroupSolverType>,
      std::unique_ptr<MomentConvergenceCheckerType>,
      std::unique_ptr<MomentCalculator>,
      const std::shared_ptr<GroupSolutionType>&,
      const UpdaterPointers& updater_ptrs,
      std::unique_ptr<MomentMapConvergenceCheckerType> moment_map_convergence_checker_ptr);
  std::unique_ptr<Initializer> BuildInitializer(
      const std::shared_ptr<formulation::updater::FixedUpdaterI>&,
      const int total_groups, const int total_angles);
  std::unique_ptr<GroupSolutionType> BuildGroupSolution(const int n_angles);
  std::unique_ptr<KEffectiveUpdaterType> BuildKEffectiveUpdater(
      const std::shared_ptr<FiniteElementType>&,
      const std::shared_ptr<CrossSectionType>&,
      const std::shared_ptr<DomainType>&);
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
      const std::shared_ptr<FissionSourceUpdater>&);
  std::unique_ptr<ParameterConvergenceCheckerType> BuildParameterConvergenceChecker(
      double max_delta, int max_iterations);
  std::shared_ptr<QuadratureSetType> BuildQuadratureSet(ParametersType);
  std::unique_ptr<SingleGroupSolverType> BuildSingleGroupSolver(
      const int max_iterations = 1000,
      const double convergence_tolerance = 1e-10);
  std::unique_ptr<SystemType> BuildSystem(const int n_groups, const int n_angles,
                                          const DomainType& domain,
                                          const std::size_t solution_size,
                                          bool is_eigenvalue_problem = true,
                                          bool need_rhs_boundary_condition = false);

 private:
  void ReportBuildingComponant(std::string componant) {
    if (!build_report_closed_) {
      data_port::StatusDataPort::Expose({"\n", utility::Color::kReset});
    }
    data_port::StatusDataPort::Expose(
        {"Building " + componant + ": ", utility::Color::kReset});
    build_report_closed_ = false;
  }

  void ReportBuildSuccess(std::string description = "") {
    data_port::StatusDataPort::Expose(
        {"Built " + description + "\n", utility::Color::kGreen});
    build_report_closed_ = true;
  }

  void Report(std::string to_report, utility::Color color) {
    data_port::StatusDataPort::Expose({to_report, color});
  }

  void ReportBuildError(std::string description = "") {
    data_port::StatusDataPort::Expose(
        {"Error: " + description + "\n", utility::Color::kRed});
    build_report_closed_ = true;
  }

  void Validate() const;

  template <typename T>
  inline std::shared_ptr<T> Shared(std::unique_ptr<T> to_convert_ptr) {
    return to_convert_ptr;
  }

  std::string ReadMappingFile(std::string filename);
  std::shared_ptr<StatusInstrument> status_instrument_ptr_{nullptr};
  std::shared_ptr<ColorStatusInstrument> color_status_instrument_ptr_{nullptr};
  std::shared_ptr<ConvergenceInstrument> convergence_status_instrument_ptr_{ nullptr };
  mutable FrameworkValidator validator_;
  const system::SystemHelper<dim> system_helper_;
  bool build_report_closed_ = true;
  std::string filename_{""};
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_
