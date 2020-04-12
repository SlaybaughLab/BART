#include "framework/builder/framework_builder.h"

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/mpi.h>
#include <sstream>

#include "utility/reporter/mpi.h"
#include "utility/reporter/colors.h"

// Convergence classes
#include "convergence/final_checker_or_n.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/reporter/mpi.h"

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

// KEffective Updater Classes
#include "calculator/cell/total_aggregated_fission_source.h"
#include "calculator/cell/integrated_fission_source.h"
#include "eigenvalue/k_effective/updater_via_fission_source.h"

// Material classes
#include "material/material_protobuf.h"

// Solver classes
#include "solver/group/single_group_solver.h"
#include "solver/gmres.h"

// Iteration classes
#include "iteration/initializer/initialize_fixed_terms_once.h"
#include "iteration/group/group_source_iteration.h"

// Quadrature classes & factories
#include "quadrature/quadrature_generator_i.h"
#include "quadrature/factory/quadrature_factories.h"
#include "quadrature/utility/quadrature_utilities.h"

// System classes
#include "system/solution/mpi_group_angular_solution.h"

namespace bart {

namespace framework {

namespace builder {

template<int dim>
void FrameworkBuilder<dim>::BuildFramework(std::string name,
                                           ParametersType& prm) {
  // Framework parameters
  int n_angles = 1; // Set to default value of 1 for scalar solve
  const int n_groups = prm.NEnergyGroups();

  *reporter_ptr_ << "Building Framework: " << Color::Green << name <<
                 Color::Reset << "\n";

  auto finite_element_ptr = Shared(BuildFiniteElement(prm));
  auto cross_sections_ptr = Shared(BuildCrossSections(prm));

  auto domain_ptr = Shared(BuildDomain(prm, finite_element_ptr,
                                       ReadMappingFile(prm.MaterialMapFilename())));
  *reporter_ptr_ << "\tSetting up domain\n";
  domain_ptr->SetUpMesh().SetUpDOF();

  std::shared_ptr<QuadratureSetType> quadrature_set_ptr = nullptr;
  UpdaterPointers updater_pointers;
  std::unique_ptr<MomentCalculatorType> moment_calculator_ptr = nullptr;

  if (prm.TransportModel() == problem::EquationType::kSelfAdjointAngularFlux) {
    quadrature_set_ptr = BuildQuadratureSet(prm);
    n_angles = quadrature_set_ptr->size();
    auto stamper_ptr = BuildStamper(domain_ptr);
    auto saaf_formulation_ptr = BuildSAAFFormulation(finite_element_ptr,
                                                     cross_sections_ptr,
                                                     quadrature_set_ptr);
    updater_pointers = BuildUpdaterPointers(
        std::move(saaf_formulation_ptr),
        std::move(stamper_ptr),
        quadrature_set_ptr);
    moment_calculator_ptr = std::move(BuildMomentCalculator(quadrature_set_ptr));

  } else if (prm.TransportModel() == problem::EquationType::kDiffusion) {
    auto diffusion_formulation_ptr = BuildDiffusionFormulation(
        finite_element_ptr,
        cross_sections_ptr);
    auto stamper_ptr = BuildStamper(domain_ptr);
    updater_pointers = BuildUpdaterPointers(
        std::move(diffusion_formulation_ptr),
        std::move(stamper_ptr));
    moment_calculator_ptr = std::move(BuildMomentCalculator());
  }

  auto initializer_ptr = BuildInitializer(
      updater_pointers.fixed_updater_ptr, n_groups, n_angles);
  auto convergence_reporter_ptr = Shared(BuildConvergenceReporter());
  auto group_solution_ptr = Shared(BuildGroupSolution(n_angles));

  auto iterative_group_solver_ptr = BuildGroupSolveIteration(
      BuildSingleGroupSolver(),
      BuildMomentConvergenceChecker(1e-10, 100),
      std::move(moment_calculator_ptr),
      group_solution_ptr,
      updater_pointers.scattering_source_updater_ptr,
      convergence_reporter_ptr);

  auto k_effective_updater = BuildKEffectiveUpdater(finite_element_ptr,
                                                    cross_sections_ptr,
                                                    domain_ptr);

  Validate();
}

template<int dim>
auto FrameworkBuilder<dim>::BuildConvergenceReporter()
-> std::unique_ptr<ReporterType> {
  reporter_ptr_->Report("\tBuilding ConvergenceReporter\n");
  using Reporter = bart::convergence::reporter::MpiNoisy;

  std::unique_ptr<ReporterType> return_ptr = nullptr;

  int this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  auto pout_ptr = std::make_unique<dealii::ConditionalOStream>(std::cout, this_process == 0);
  return_ptr = std::make_unique<Reporter>(std::move(pout_ptr));

  return std::move(return_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildCrossSections(
    const problem::ParametersI& problem_parameters)
    -> std::unique_ptr<CrossSectionType> {
  reporter_ptr_->Report("\tBuilding Cross-sections: ");
  std::unique_ptr<CrossSectionType> return_ptr = nullptr;
  // Default implementation using protocol buffers
  try {
    MaterialProtobuf materials(problem_parameters.MaterialFilenames(),
                               problem_parameters.IsEigenvalueProblem(),
                               problem_parameters.DoNDA(),
                               problem_parameters.NEnergyGroups(),
                               problem_parameters.NumberOfMaterials());
    return_ptr = std::move(std::make_unique<CrossSectionType>(materials));
    reporter_ptr_->Report("Built (default) Cross-sections using protobuf\n",
                          utility::reporter::Color::Green);
  } catch (...) {
    reporter_ptr_->Report("Error building (default) Cross-sections using protobuf\n",
                          utility::reporter::Color::Red);
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
    std::unique_ptr<StamperType> stamper_ptr)
-> UpdaterPointers {
  ReportBuildingComponant("Building Diffusion Formulation updater");
  UpdaterPointers return_struct;

  using ReturnType = formulation::updater::DiffusionUpdater<dim>;

  auto diffusion_updater_ptr = std::make_shared<ReturnType>(
      std::move(formulation_ptr),
      std::move(stamper_ptr));
  return_struct.fixed_updater_ptr = diffusion_updater_ptr;
  return_struct.scattering_source_updater_ptr = diffusion_updater_ptr;
  return_struct.fission_source_updater_ptr = diffusion_updater_ptr;
  ReportBuildSuccess("");
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
  return_struct.fixed_updater_ptr = saaf_updater_ptr;
  return_struct.scattering_source_updater_ptr = saaf_updater_ptr;
  return_struct.fission_source_updater_ptr = saaf_updater_ptr;
  ReportBuildSuccess("");
  return return_struct;
}

template <int dim>
auto FrameworkBuilder<dim>::BuildGroupSolveIteration(
    std::unique_ptr<SingleGroupSolverType> single_group_solver_ptr,
    std::unique_ptr<MomentConvergenceCheckerType> moment_convergence_checker_ptr,
    std::unique_ptr<MomentCalculatorType> moment_calculator_ptr,
    const std::shared_ptr<GroupSolutionType>& group_solution_ptr,
    const std::shared_ptr<ScatteringSourceUpdaterType>& scattering_source_updater_ptr,
    const std::shared_ptr<ReporterType>& convergence_report_ptr)
    -> std::unique_ptr<GroupSolveIterationType> {
  std::unique_ptr<GroupSolveIterationType> return_ptr = nullptr;

  ReportBuildingComponant("Iterative group solver");

  return_ptr = std::move(
      std::make_unique<iteration::group::GroupSourceIteration<dim>>(
          std::move(single_group_solver_ptr),
          std::move(moment_convergence_checker_ptr),
          std::move(moment_calculator_ptr),
          group_solution_ptr,
          scattering_source_updater_ptr,
          convergence_report_ptr)
      );
  has_scattering_source_update_ = true;
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
  reporter_ptr_->Report("\tBuilding Initializer\n");
  std::unique_ptr<InitializerType> return_ptr = nullptr;

  using InitializeOnceType = iteration::initializer::InitializeFixedTermsOnce;

  return_ptr = std::move(std::make_unique<InitializeOnceType>(
      updater_ptr, total_groups, total_angles));

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

  ReportBuildSuccess("");
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
  reporter_ptr_->Report("\tBuilding Moment Calculator: ");
  std::unique_ptr<MomentCalculatorType> return_ptr = nullptr;

  try {
    return_ptr = std::move(quadrature::factory::MakeMomentCalculator<dim>(
        implementation, quadrature_set_ptr));

    if (implementation == MomentCalculatorImpl::kScalarMoment) {
      reporter_ptr_->Report("Built (default) calculator for scalar solve\n",
                            Color::Green);
    } else if (implementation == MomentCalculatorImpl::kZerothMomentOnly) {
      reporter_ptr_->Report("Built (default) calculator for 0th moment only\n",
                            Color::Green);
    } else {
      AssertThrow(false,
                  dealii::ExcMessage("Unsupported implementation of moment "
                                     "calculator specified in call to "
                                     "BuildMomentCalculator"))
    }
  } catch (...) {
    reporter_ptr_->Report("Error building calculator for scalar solve \n",
                          Color::Red);
    throw;
  }

  return return_ptr;
}



template<int dim>
auto FrameworkBuilder<dim>::BuildMomentConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentConvergenceCheckerType>{
  //TODO(Josh): Add option for using other than L1Norm
  reporter_ptr_->Report("\tBuilding Moment Convergence Checker\n");
  using CheckerType = convergence::moments::SingleMomentCheckerL1Norm;
  using FinalCheckerType = convergence::FinalCheckerOrN<
      system::moments::MomentVector,
      convergence::moments::SingleMomentCheckerI>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));

  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
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
  ReportBuildSuccess("");
  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
}

template<int dim>
auto FrameworkBuilder<dim>::BuildQuadratureSet(ParametersType problem_parameters)
-> std::shared_ptr<QuadratureSetType> {
  reporter_ptr_->Report("\tBuilding quadrature set\n");
  using QuadratureGeneratorType = quadrature::QuadratureGeneratorI<dim>;

  std::shared_ptr<QuadratureSetType> return_ptr = nullptr;
  std::shared_ptr<QuadratureGeneratorType > quadrature_generator_ptr = nullptr;

  const int order_value = problem_parameters.AngularQuadOrder();
  switch (problem_parameters.AngularQuad()) {
    default: {
      if (dim == 3) {
        quadrature_generator_ptr =
            quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
                quadrature::Order(order_value),
                quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
      } else {
        AssertThrow(false,
                    dealii::ExcMessage("No supported quadratures for this dimension "
                                       "and transport model"))
      }
    }
  }

  return_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

  auto quadrature_points = quadrature::utility::GenerateAllPositiveX<dim>(
      quadrature_generator_ptr->GenerateSet());

  quadrature::factory::FillQuadratureSet<dim>(return_ptr.get(),
                                              quadrature_points);

  return std::move(return_ptr);
}

template <int dim>
auto FrameworkBuilder<dim>::BuildSAAFFormulation(
    const std::shared_ptr<FiniteElementType>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const std::shared_ptr<QuadratureSetType>& quadrature_set_ptr,
    const formulation::SAAFFormulationImpl implementation)
-> std::unique_ptr<SAAFFormulationType> {
  reporter_ptr_->Report("\tBuilding SAAF Formulation\n");
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
  reporter_ptr_->Report("\tBuilding SingleGroupSolver: ");
  std::unique_ptr<SingleGroupSolverType> return_ptr = nullptr;

  auto linear_solver_ptr = std::make_unique<solver::GMRES>(max_iterations,
                                                           convergence_tolerance);
  reporter_ptr_->Report("Built GMRES: tol = " + std::to_string(convergence_tolerance)
                            + "iter_max = " + std::to_string(max_iterations) + "\n",
                            Color::Green);
  return_ptr = std::move(std::make_unique<solver::group::SingleGroupSolver>(
          std::move(linear_solver_ptr)));

  return return_ptr;
}

template<int dim>
auto FrameworkBuilder<dim>::BuildStamper(
    const std::shared_ptr<DomainType>& domain_ptr)
-> std::unique_ptr<StamperType> {
  reporter_ptr_->Report("\tBuilding Stamper\n");
  std::unique_ptr<StamperType> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<formulation::Stamper<dim>>(domain_ptr));

  return return_ptr;
}
template<int dim>
std::string FrameworkBuilder<dim>::ReadMappingFile(std::string filename) {
  reporter_ptr_->Report("\tReading mapping file: ");

  std::ifstream mapping_file(filename);
  if (mapping_file.is_open()) {
    reporter_ptr_->Report(filename + '\n', Color::Green);

    return std::string(
        (std::istreambuf_iterator<char>(mapping_file)),
        std::istreambuf_iterator<char>());
  } else {
    reporter_ptr_->Report("Error reading " + filename + "\n", Color::Red);
    AssertThrow(false,
                dealii::ExcMessage("Failed to open material mapping file"))
  }
}

template<int dim>
void FrameworkBuilder<dim>::Validate() const {
  bool issue = false;
  reporter_ptr_->Report("Validating framework\n");
  reporter_ptr_->Report("\tHas scattering source update: ");
  if (has_scattering_source_update_) {
    reporter_ptr_->Report("True\n", Color::Green);
  } else {
    reporter_ptr_->Report("False\n", Color::Red);
    issue = true;
  }
  reporter_ptr_->Report("\tHas fission source update: ");
  if (has_fission_source_update_) {
    reporter_ptr_->Report("True\n", Color::Green);
  } else {
    reporter_ptr_->Report("False\n", Color::Red);
    issue = true;
  }
  if (issue) {
    reporter_ptr_->Report("Warning: one or more issues identified during "
                          "framework validation\n", Color::Yellow);
  }
}

template class FrameworkBuilder<1>;
template class FrameworkBuilder<2>;
template class FrameworkBuilder<3>;

} // namespace builder

} // namespace framework

} // namespace bart
