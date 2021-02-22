#ifndef BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_
#define BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_

#include "framework/builder/framework_builder_i.hpp"

#include <deal.II/base/conditional_ostream.h>

#include "framework/builder/framework_validator.hpp"
#include "instrumentation/port.hpp"
#include "system/system_helper.hpp"
#include "utility/has_description.h"


namespace bart::framework::builder {

namespace data_port {
struct BuilderStatus;
using StatusDataPort = instrumentation::Port<std::pair<std::string, utility::Color>, BuilderStatus>;
}

template <int dim>
class FrameworkBuilder : public data_port::StatusDataPort, public FrameworkBuilderI<dim> {
 public:
  // New using types from refactor
  using typename FrameworkBuilderI<dim>::AngularFluxIntegrator;
  using typename FrameworkBuilderI<dim>::CrossSections;
  using typename FrameworkBuilderI<dim>::DiffusionFormulation;
  using typename FrameworkBuilderI<dim>::DriftDiffusionFormulation;
  using typename FrameworkBuilderI<dim>::Domain;
  using typename FrameworkBuilderI<dim>::FiniteElement;
  using typename FrameworkBuilderI<dim>::FrameworkI;
  using typename FrameworkBuilderI<dim>::GroupSolution;
  using typename FrameworkBuilderI<dim>::GroupSolveIteration;
  using typename FrameworkBuilderI<dim>::Initializer;
  using typename FrameworkBuilderI<dim>::KEffectiveUpdater;
  using typename FrameworkBuilderI<dim>::MomentCalculator;
  using typename FrameworkBuilderI<dim>::MomentConvergenceChecker;
  using typename FrameworkBuilderI<dim>::MomentMapConvergenceChecker;
  using typename FrameworkBuilderI<dim>::OuterIteration;
  using typename FrameworkBuilderI<dim>::ParameterConvergenceChecker;
  using typename FrameworkBuilderI<dim>::QuadratureSet;
  using typename FrameworkBuilderI<dim>::Stamper;
  using typename FrameworkBuilderI<dim>::Subroutine;
  using typename FrameworkBuilderI<dim>::SubroutineName;
  using typename FrameworkBuilderI<dim>::SAAFFormulation;
  using typename FrameworkBuilderI<dim>::SphericalHarmonicMoments;
  using typename FrameworkBuilderI<dim>::SingleGroupSolver;
  using typename FrameworkBuilderI<dim>::System;
  using typename FrameworkBuilderI<dim>::Validator;

  using typename FrameworkBuilderI<dim>::UpdaterPointers;
  using typename FrameworkBuilderI<dim>::BoundaryConditionsUpdater;
  using typename FrameworkBuilderI<dim>::FissionSourceUpdater;
  using typename FrameworkBuilderI<dim>::FixedTermUpdater;
  using typename FrameworkBuilderI<dim>::ScatteringSourceUpdater;

  using typename FrameworkBuilderI<dim>::DiffusionFormulationImpl;
  using typename FrameworkBuilderI<dim>::MomentCalculatorImpl;
  using typename FrameworkBuilderI<dim>::InitializerName;

  using typename FrameworkBuilderI<dim>::ColorStatusPair;
  using typename FrameworkBuilderI<dim>::ColorStatusInstrument;
  using typename FrameworkBuilderI<dim>::ConvergenceInstrument;
  using typename FrameworkBuilderI<dim>::StatusInstrument;

  using typename FrameworkBuilderI<dim>::AngularFluxStorage;

  FrameworkBuilder(std::unique_ptr<Validator> validator_ptr);
  ~FrameworkBuilder() = default;

  [[nodiscard]] auto BuildAngularFluxIntegrator(
      const std::shared_ptr<QuadratureSet>) -> std::unique_ptr<AngularFluxIntegrator> override;

  [[nodiscard]] auto BuildDiffusionFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::cross_sections::CrossSections>&,
      const DiffusionFormulationImpl implementation = DiffusionFormulationImpl::kDefault)
  -> std::unique_ptr<DiffusionFormulation> override;

  [[nodiscard]] auto BuildDriftDiffusionFormulation(
      const std::shared_ptr<AngularFluxIntegrator>&,
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::cross_sections::CrossSections>&) -> std::unique_ptr<DriftDiffusionFormulation> override;

  [[nodiscard]] auto BuildDomain(const FrameworkParameters::DomainSize,
                                 const FrameworkParameters::NumberOfCells,
                                 const std::shared_ptr<FiniteElement>&,
                                 const std::string material_mapping) -> std::unique_ptr<Domain> override;
  [[nodiscard]] auto BuildFiniteElement(
      const problem::CellFiniteElementType finite_element_type,
      const problem::DiscretizationType discretization_type,
      const FrameworkParameters::PolynomialDegree polynomial_degree) -> std::unique_ptr<FiniteElement> override;
  [[nodiscard]] auto BuildGroupSolution(const int n_angles) -> std::unique_ptr<GroupSolution> override;
  [[nodiscard]] auto BuildGroupSolveIteration(
      std::unique_ptr<SingleGroupSolver>,
      std::unique_ptr<MomentConvergenceChecker>,
      std::unique_ptr<MomentCalculator>,
      const std::shared_ptr<GroupSolution>&,
      const UpdaterPointers& updater_ptrs,
      std::unique_ptr<MomentMapConvergenceChecker>) -> std::unique_ptr<GroupSolveIteration> override;
  [[nodiscard]] auto BuildInitializer(const std::shared_ptr<FixedTermUpdater>&,
                                      const int total_groups,
                                      const int total_angles) -> std::unique_ptr<Initializer> override;
  [[nodiscard]] auto BuildInitializer(const std::shared_ptr<FixedTermUpdater>&,
                                      const int total_groups,
                                      const int total_angles,
                                      const InitializerName) -> std::unique_ptr<Initializer> override;
  [[nodiscard]] auto BuildKEffectiveUpdater() -> std::unique_ptr<KEffectiveUpdater> override;
  [[nodiscard]] auto BuildKEffectiveUpdater(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<CrossSections>&,
      const std::shared_ptr<Domain>&) -> std::unique_ptr<KEffectiveUpdater> override;
  [[nodiscard]] auto BuildMomentCalculator(
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kScalarMoment)
  -> std::unique_ptr<MomentCalculator> override;
  [[nodiscard]] auto BuildMomentCalculator(
      std::shared_ptr<QuadratureSet>,
      MomentCalculatorImpl implementation = MomentCalculatorImpl::kZerothMomentOnly)
  -> std::unique_ptr<MomentCalculator> override;
  [[nodiscard]] auto  BuildMomentConvergenceChecker(
      double max_delta,
      int max_iterations) -> std::unique_ptr<MomentConvergenceChecker> override;
  [[nodiscard]] auto BuildMomentMapConvergenceChecker(
      double max_delta,
      int max_iterations) -> std::unique_ptr<MomentMapConvergenceChecker> override;
  [[nodiscard]] auto BuildOuterIteration(
      std::unique_ptr<GroupSolveIteration>,
      std::unique_ptr<ParameterConvergenceChecker>,
      const std::string& output_filename_base) -> std::unique_ptr<OuterIteration> override;
  [[nodiscard]] auto BuildOuterIteration(
      std::unique_ptr<GroupSolveIteration>,
      std::unique_ptr<ParameterConvergenceChecker>,
      std::unique_ptr<KEffectiveUpdater>,
      const std::shared_ptr<FissionSourceUpdater>&,
      const std::string& output_filename_base) -> std::unique_ptr<OuterIteration> override;
  [[nodiscard]] auto BuildParameterConvergenceChecker(
      double max_delta,
      int max_iterations) -> std::unique_ptr<ParameterConvergenceChecker> override;
  [[nodiscard]] auto BuildQuadratureSet(
      const problem::AngularQuadType,
      const FrameworkParameters::AngularQuadratureOrder) -> std::shared_ptr<QuadratureSet> override;
  [[nodiscard]] auto BuildSAAFFormulation(
      const std::shared_ptr<FiniteElement>&,
      const std::shared_ptr<data::cross_sections::CrossSections>&,
      const std::shared_ptr<QuadratureSet>&,
      const formulation::SAAFFormulationImpl implementation = formulation::SAAFFormulationImpl::kDefault)
  -> std::unique_ptr<SAAFFormulation> override;
  [[nodiscard]] auto BuildSingleGroupSolver(
      const int max_iterations,
      const double convergence_tolerance) -> std::unique_ptr<SingleGroupSolver> override;
  [[nodiscard]] auto BuildStamper(const std::shared_ptr<Domain>&) -> std::unique_ptr<Stamper> override;
  [[nodiscard]] auto BuildSubroutine(std::unique_ptr<FrameworkI>,
                                     const SubroutineName) -> std::unique_ptr<Subroutine> override;
  [[nodiscard]] auto BuildSystem(const int n_groups,
                                 const int n_angles,
                                 const Domain& domain,
                                 const std::size_t solution_size,
                                 bool is_eigenvalue_problem,
                                 bool need_rhs_boundary_condition) -> std::unique_ptr<System> override;

  [[nodiscard]] auto BuildUpdaterPointers(std::unique_ptr<DiffusionFormulation>,
                                    std::unique_ptr<DriftDiffusionFormulation>,
                                    std::shared_ptr<Stamper>,
                                    std::shared_ptr<AngularFluxIntegrator>,
                                    std::shared_ptr<SphericalHarmonicMoments>,
                                    AngularFluxStorage&,
                                    const std::map<problem::Boundary, bool>&) -> UpdaterPointers override;
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

  // Instrumentation
  auto set_color_status_instrument_ptr(
      const std::shared_ptr<ColorStatusInstrument>& to_set) -> FrameworkBuilderI<dim>& override {
    AssertThrow(to_set != nullptr, dealii::ExcMessage("Error setting instrument, pointer is null"))
    color_status_instrument_ptr_ = to_set;
    return *this;
  };
  auto set_convergence_status_instrument_ptr(
      const std::shared_ptr<ConvergenceInstrument>& to_set) -> FrameworkBuilderI<dim>& override {
    AssertThrow(to_set != nullptr, dealii::ExcMessage("Error setting instrument, pointer is null"))
    convergence_status_instrument_ptr_ = to_set;
    return *this;
  };
  auto set_status_instrument_ptr(const std::shared_ptr<StatusInstrument>& to_set) -> FrameworkBuilderI<dim>& override{
    AssertThrow(to_set != nullptr, dealii::ExcMessage("Error setting instrument, pointer is null"))
    status_instrument_ptr_ = to_set;
    return *this;
  };

  auto color_status_instrument_ptr() const -> std::shared_ptr<ColorStatusInstrument> override {
    return color_status_instrument_ptr_; };
  auto convergence_status_instrument_ptr() const -> std::shared_ptr<ConvergenceInstrument> override {
    return convergence_status_instrument_ptr_; };
  auto status_instrument_ptr() const -> std::shared_ptr<StatusInstrument> override {
    return status_instrument_ptr_; };

  auto validator_ptr() -> Validator* override { return validator_ptr_.get(); };

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
  std::unique_ptr<Validator> validator_ptr_{ nullptr };
  std::shared_ptr<StatusInstrument> status_instrument_ptr_{nullptr};
  std::shared_ptr<ColorStatusInstrument> color_status_instrument_ptr_{nullptr};
  std::shared_ptr<ConvergenceInstrument> convergence_status_instrument_ptr_{ nullptr };
  const system::SystemHelper<dim> system_helper_;
  bool build_report_closed_ = true;
};

} // namespace bart::framework::builder

#endif //BART_SRC_FRAMEWORK_BUILDER_FRAMEWORK_BUILDER_HPP_
