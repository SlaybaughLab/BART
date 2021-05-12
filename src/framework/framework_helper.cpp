#include "framework/framework_helper.hpp"

#include "solver/eigenvalue/krylov_schur_eigenvalue_solver.hpp"
#include "acceleration/two_grid/spectral_shape/domain_spectral_shapes.hpp"
#include "acceleration/two_grid/spectral_shape/material_spectral_shapes.hpp"
#include "acceleration/two_grid/spectral_shape/spectral_shape.hpp"
#include "acceleration/two_grid/flux_corrector.hpp"
#include "calculator/residual/cell_isotropic_residual.hpp"
#include "calculator/residual/domain_isotropic_residual.hpp"
#include "formulation/scalar/two_grid_diffusion.hpp"

#include "framework/builder/framework_builder.hpp"
#include "framework/framework.hpp"
#include "framework/builder/framework_validator.hpp"
#include "instrumentation/builder/instrument_builder.hpp"
#include "instrumentation/converter/convert_to_string/convergence_to_string.h"
#include "data/material/material_protobuf.hpp"
#include "data/cross_sections/collapsed_one_group_cross_sections.hpp"
#include "iteration/outer/outer_iteration.hpp"
#include "results/output_dealii_vtu.h"
#include "system/system_helper.hpp"
#include "system/solution/solution_types.h"

// to be removed
#include "instrumentation/basic_instrument.h"
#include "instrumentation/outstream/vector_to_vtu.hpp"
#include "iteration/outer/outer_fixed_source_iteration.hpp"
#include "iteration/group/group_source_iteration.hpp"
#include "formulation/updater/fixed_updater.hpp"
#include "iteration/subroutine/two_grid_acceleration.hpp"


#include <fstream>

#include <fmt/color.h>
#include <system/system_helper.hpp>

namespace bart::framework {

namespace  {
using InstrumentBuilder = instrumentation::builder::InstrumentBuilder;
using InstrumentName = instrumentation::builder::InstrumentName;
template <typename T> inline std::shared_ptr<T> Shared(std::unique_ptr<T> to_convert_ptr) { return to_convert_ptr; }
} // namespace

template<int dim>
FrameworkHelper<dim>::FrameworkHelper(const std::shared_ptr<SystemHelper>& system_helper_ptr)
    : system_helper_ptr_(system_helper_ptr) {
  AssertThrow(system_helper_ptr_ != nullptr, dealii::ExcMessage("Error in constructor of framework helper, system"
                                                                "helper pointer passed is null"))
}

template<int dim>
auto FrameworkHelper<dim>::ToFrameworkParameters(
    const problem::ParametersI &problem_parameters) -> framework::FrameworkParameters {
  using Boundary = problem::Boundary;

  framework::FrameworkParameters return_parameters {
    .output_filename_base{ problem_parameters.OutputFilenameBase()},
    .neutron_energy_groups{ problem_parameters.NEnergyGroups() },
    .equation_type{ problem_parameters.TransportModel() },
    .k_effective_updater{ problem_parameters.K_EffectiveUpdaterType() },
    .group_solver_type{ problem_parameters.InGroupSolver() },
    .angular_quadrature_type{ problem_parameters.AngularQuad() },
    .angular_quadrature_order{ quadrature::Order(problem_parameters.AngularQuadOrder()) },
    .spatial_dimension{ framework::FrameworkParameters::SpatialDimension(problem_parameters.SpatialDimension()) },
    .domain_size{ framework::FrameworkParameters::DomainSize(problem_parameters.SpatialMax()) },
    .number_of_cells{ framework::FrameworkParameters::NumberOfCells(problem_parameters.NCells()) },
    .uniform_refinements{ problem_parameters.UniformRefinements() },
    .discretization_type{ problem_parameters.Discretization() },
    .polynomial_degree{ framework::FrameworkParameters::PolynomialDegree(problem_parameters.FEPolynomialDegree()) },
    .use_nda_{ problem_parameters.DoNDA() },
    .use_two_grid_{ problem_parameters.UseTwoGridAcceleration() },
    .output_aggregated_source_data{ problem_parameters.OutputAggregatedSourceData() },
    .output_scalar_flux_as_vtu{ problem_parameters.OutputScalarFluxAsVTU() },
    .output_fission_source_as_vtu{ problem_parameters.OutputFissionSourceAsVTU() },
    .output_scattering_source_as_vtu{ problem_parameters.OutputScatteringSourceAsVTU() },
    .output_inner_iterations_to_file{ problem_parameters.OutputInnerIterationsToFile() }
  };

  std::set<Boundary> reflective_boundaries;
  for (const auto& [boundary, is_reflective] : problem_parameters.ReflectiveBoundary()){
    if (is_reflective)
      reflective_boundaries.insert(boundary);
  }
  return_parameters.reflective_boundaries = reflective_boundaries;

  // Open material mapping file and read
  std::ifstream mapping_file(problem_parameters.MaterialMapFilename());
  if (mapping_file.is_open()) {
    return_parameters.material_mapping = std::string((std::istreambuf_iterator<char>(mapping_file)),
                                                     std::istreambuf_iterator<char>());
  } else {
    AssertThrow(false, dealii::ExcMessage("Failed to open material mapping file."))
  }

  const auto eigen_solver_type{ problem_parameters.EigenSolver() };
  const bool is_eigenvalue_solve{ eigen_solver_type != problem::EigenSolverType::kNone };

  if (is_eigenvalue_solve) {
    return_parameters.eigen_solver_type = eigen_solver_type;
  }

  data::material::MaterialProtobuf materials(problem_parameters.MaterialFilenames(),
                                       is_eigenvalue_solve,
                                       false,
                                       return_parameters.neutron_energy_groups,
                                       problem_parameters.NumberOfMaterials());

  return_parameters.cross_sections_ = std::make_shared<data::cross_sections::MaterialCrossSections>(materials);

  return return_parameters;
}

template<int dim>
auto FrameworkHelper<dim>::BuildFramework(
    builder::FrameworkBuilderI<dim>& builder,
    framework::FrameworkParameters& parameters) -> std::unique_ptr<framework::FrameworkI> {
  using FrameworkPart = framework::builder::FrameworkPart;
  using MomentCalculator = typename builder::FrameworkBuilderI<dim>::MomentCalculator;
  using OuterIteration = typename builder::FrameworkBuilderI<dim>::OuterIteration;
  using MomentCalculatorImpl = typename builder::FrameworkBuilderI<dim>::MomentCalculatorImpl;
  using UpdaterPointers = typename builder::FrameworkBuilderI<dim>::UpdaterPointers;
  using QuadratureSet = typename builder::FrameworkBuilderI<dim>::QuadratureSet;
  using Validator = framework::builder::FrameworkValidator;

  using ColorStringPair = std::pair<std::string, utility::Color>;

  // Check that if Drift-Diffusion formujlation is specified that we have the required dependencies provided
  if (parameters.equation_type == problem::EquationType::kDriftDiffusion) {
    std::string nda_dependency_error{ "Error building framework for Drift Diffusion: Required dependency missing: "};
    AssertThrow(parameters.nda_data_.angular_flux_integrator_ptr_ != nullptr,
        dealii::ExcMessage(nda_dependency_error + "angular flux integrator ptr is null"))
    AssertThrow(parameters.nda_data_.higher_order_moments_ptr_ != nullptr,
                dealii::ExcMessage(nda_dependency_error + "higher order moments ptr is null"))
  }



  // Build instruments to be used
  auto color_string_instrument_ptr = Shared(
      InstrumentBuilder::BuildInstrument<ColorStringPair>(InstrumentName::kColorStatusToConditionalOstream));
  auto convergence_status_instrument_ptr = Shared(
      InstrumentBuilder::BuildInstrument<convergence::Status>(InstrumentName::kConvergenceStatusToConditionalOstream));
  auto string_instrument_ptr = Shared(
      InstrumentBuilder::BuildInstrument<std::string>(InstrumentName::kStringToConditionalOstream));

  builder.set_color_status_instrument_ptr(color_string_instrument_ptr);
  builder.set_convergence_status_instrument_ptr(convergence_status_instrument_ptr);
  builder.set_status_instrument_ptr(string_instrument_ptr);

  try {
    auto& dynamic_framework_builder = dynamic_cast<framework::builder::FrameworkBuilder<dim>&>(builder);
    instrumentation::GetPort<framework::builder::data_port::StatusDataPort>(dynamic_framework_builder)
        .AddInstrument(color_string_instrument_ptr);
  } catch (std::bad_cast&) {}

  auto& validator = *builder.validator_ptr();
  validator.Parse(parameters);

  const int n_groups{ parameters.neutron_energy_groups };
  const bool need_angular_solution_storage{ validator.NeededParts().contains(FrameworkPart::AngularSolutionStorage) };
  const bool has_reflective_boundaries { !parameters.reflective_boundaries.empty() };
  const std::string output_filename_base { parameters.output_filename_base };

  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "Building framework: {}\n", parameters.name);

  auto finite_element_ptr = Shared(builder.BuildFiniteElement(parameters.cell_finite_element_type,
                                                              parameters.discretization_type,
                                                              parameters.polynomial_degree));
  auto domain_ptr = Shared(builder.BuildDomain(parameters.domain_size, parameters.number_of_cells,
                                               finite_element_ptr, parameters.material_mapping));

  fmt::print("Setting up domain...\n");
  domain_ptr->SetUpMesh(parameters.uniform_refinements);
  domain_ptr->SetUpDOF();

  // These objects will be set up differently depending on the implementation
  std::shared_ptr<QuadratureSet> quadrature_set_ptr{ nullptr };
  UpdaterPointers updater_pointers;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr{ nullptr };


  // Set up for Angular/Scalar solve
  using EquationType = problem::EquationType;
  std::set<EquationType> scalar_types{EquationType::kDiffusion, EquationType::kDriftDiffusion};
  std::set<EquationType> angular_types{EquationType::kSelfAdjointAngularFlux};

  if (angular_types.contains(parameters.equation_type)) {
    // Angular solve
    AssertThrow(parameters.angular_quadrature_order.has_value(),
                dealii::ExcMessage("Error building framework, equation type requires quadrature but order is null"))
    quadrature_set_ptr = builder.BuildQuadratureSet(parameters.angular_quadrature_type,
                                                    parameters.angular_quadrature_order.value());
    moment_calculator_ptr = builder.BuildMomentCalculator(quadrature_set_ptr, MomentCalculatorImpl::kZerothMomentOnly);
  } else {
    // Scalar solve
    moment_calculator_ptr = builder.BuildMomentCalculator(MomentCalculatorImpl::kScalarMoment);
  }

  const int n_angles { quadrature_set_ptr == nullptr ? 1 : static_cast<int>(quadrature_set_ptr->size()) };

  // Set up angular solutions if needed
  system::solution::EnergyGroupToAngularSolutionPtrMap angular_solutions_;
  if (need_angular_solution_storage)
    system_helper_ptr_->SetUpEnergyGroupToAngularSolutionPtrMap(angular_solutions_, n_groups, n_angles);

  //TODO: Add overload that makes this unecessary
  std::map<problem::Boundary, bool> reflective_boundaries {
      {problem::Boundary::kXMin, false}, {problem::Boundary::kXMax, false},
      {problem::Boundary::kYMin, false}, {problem::Boundary::kYMax, false},
      {problem::Boundary::kZMin, false}, {problem::Boundary::kZMax, false},
  };
  for (auto& boundary : parameters.reflective_boundaries) {
    reflective_boundaries.at(boundary) = true;
  }

  if (parameters.equation_type == problem::EquationType::kTwoGridDiffusion) {
    AssertThrow(parameters.two_grid_data_.one_group_cross_sections_ptr_ != nullptr,
                dealii::ExcMessage("Not one-group cross_sections"));

    auto two_grid_diffusion_formulation = std::make_unique<formulation::scalar::TwoGridDiffusion<dim>>(
        finite_element_ptr, parameters.cross_sections_.value(),
        parameters.two_grid_data_.one_group_cross_sections_ptr_);
    two_grid_diffusion_formulation->Precalculate(domain_ptr->Cells().at(0));
    updater_pointers = builder.BuildUpdaterPointers(std::move(two_grid_diffusion_formulation),
                                                    builder.BuildStamper(domain_ptr),
                                                    reflective_boundaries);
    std::cout << "Build two-grid diffusion formulation\n";
  } else if (parameters.equation_type == problem::EquationType::kSelfAdjointAngularFlux) {
    auto saaf_formulation_ptr = builder.BuildSAAFFormulation(finite_element_ptr,
                                                             parameters.cross_sections_.value(),
                                                             quadrature_set_ptr,
                                                             formulation::SAAFFormulationImpl::kDefault);
    saaf_formulation_ptr->Initialize(domain_ptr->Cells().at(0));
    if (has_reflective_boundaries) {
      updater_pointers = builder.BuildUpdaterPointers(std::move(saaf_formulation_ptr),
                                                      builder.BuildStamper(domain_ptr),
                                                      quadrature_set_ptr,
                                                      reflective_boundaries,
                                                      angular_solutions_);
    } else {
      updater_pointers = builder.BuildUpdaterPointers(std::move(saaf_formulation_ptr),
                                                      builder.BuildStamper(domain_ptr),
                                                      quadrature_set_ptr);
    }
  } else if (parameters.equation_type == EquationType::kDiffusion ||
      parameters.equation_type == EquationType::kDriftDiffusion) {
    auto diffusion_formulation_ptr = builder.BuildDiffusionFormulation(finite_element_ptr,
                                                                       parameters.cross_sections_.value(),
                                                                       formulation::DiffusionFormulationImpl::kDefault);
    diffusion_formulation_ptr->Precalculate(domain_ptr->Cells().at(0));
    if (parameters.equation_type == EquationType::kDriftDiffusion) {
      auto drift_diffusion_formulation_ptr = builder.BuildDriftDiffusionFormulation(
          parameters.nda_data_.angular_flux_integrator_ptr_,
          finite_element_ptr,
          parameters.cross_sections_.value());
      updater_pointers = builder.BuildUpdaterPointers(
          std::move(diffusion_formulation_ptr),
          std::move(drift_diffusion_formulation_ptr),
          builder.BuildStamper(domain_ptr),
          parameters.nda_data_.angular_flux_integrator_ptr_,
          parameters.nda_data_.higher_order_moments_ptr_,
          parameters.nda_data_.higher_order_angular_flux_,
          reflective_boundaries);
    } else {
      updater_pointers = builder.BuildUpdaterPointers(std::move(diffusion_formulation_ptr),
                                                      builder.BuildStamper(domain_ptr),
                                                      reflective_boundaries);
    }
  }

  if (parameters.output_aggregated_source_data) {
    auto make_source_instrument = [](const std::string filename) {
      return Shared(InstrumentBuilder::BuildInstrument<double>(
          InstrumentName::kDoubleToFile, filename));
    };
    if (updater_pointers.fission_source_updater_ptr != nullptr) {
      instrumentation::GetPort<formulation::updater::data_port::AggregatedFissionSourceValue>(
          *updater_pointers.fission_source_updater_ptr.get())
          .AddInstrument(make_source_instrument(parameters.output_filename_base + "_fission_source.csv"));
    }
    if (updater_pointers.scattering_source_updater_ptr != nullptr) {
      instrumentation::GetPort<formulation::updater::data_port::AggregatedScatteringSourceValue>(
          *updater_pointers.scattering_source_updater_ptr.get())
          .AddInstrument(make_source_instrument(parameters.output_filename_base + "_scattering_source.csv"));
    }
    if (updater_pointers.boundary_conditions_updater_ptr != nullptr) {
      instrumentation::GetPort<formulation::updater::data_port::AggregatedBoundaryConditionValue>(
          *updater_pointers.boundary_conditions_updater_ptr.get())
          .AddInstrument(make_source_instrument(parameters.output_filename_base + "_boundary_source.csv"));
    }
  }

  auto initializer_ptr = builder.BuildInitializer(updater_pointers.fixed_updater_ptr,
                                                  parameters.neutron_energy_groups,
                                                  n_angles,
                                                  iteration::initializer::InitializerName::kInitializeFixedTermsAndResetMoments);

  auto group_solution_ptr = Shared(builder.BuildGroupSolution(n_angles));
  system_helper_ptr_->SetUpMPIAngularSolution(*group_solution_ptr, *domain_ptr, 1.0);

  auto group_iteration_ptr = builder.BuildGroupSolveIteration(
      builder.BuildSingleGroupSolver(10000, 1e-10),
      builder.BuildMomentConvergenceChecker(1e-6, 1000),
      std::move(moment_calculator_ptr),
      group_solution_ptr,
      updater_pointers,
      builder.BuildMomentMapConvergenceChecker(1e-6, 1000));

  if (parameters.output_inner_iterations_to_file) {
    try {
      instrumentation::GetPort<iteration::group::data_ports::NumberOfIterationsPort>(*group_iteration_ptr)
          .AddInstrument(Shared(InstrumentBuilder::BuildInstrument<double>(InstrumentName::kDoubleToFile,
                                                                           parameters.output_filename_base
                                                                               + "_inner_iterations.csv")));
    } catch (std::bad_cast&) {
      std::cout << "Warning: Output Inner Iterations to file was selected but constructed group iteration class does "
                   "not support required instrumentation.\n";
    }
  }

  if (need_angular_solution_storage) {
    group_iteration_ptr->UpdateThisAngularSolutionMap(angular_solutions_);
    validator.AddPart(FrameworkPart::AngularSolutionStorage);
  }

  auto system_ptr = builder.BuildSystem(parameters.neutron_energy_groups,
                                        n_angles,
                                        *domain_ptr,
                                        group_solution_ptr->GetSolution(0).size(),
                                        parameters.eigen_solver_type.has_value(),
                                        need_angular_solution_storage);

  std::unique_ptr<iteration::subroutine::SubroutineI> group_post_processing_subroutine{ nullptr };
  if (parameters.use_two_grid_) {
    std::cout << "Setting up two_grid ===============================================================================\n";
    acceleration::two_grid::spectral_shape::MaterialSpectralShapes material_spectral_shape_calculator(
        std::make_unique<acceleration::two_grid::spectral_shape::SpectralShape>(
            std::make_unique<solver::eigenvalue::KrylovSchurEigenvalueSolver>()));
    std::cout << "calculate material Spectral Shape" << std::endl;
    material_spectral_shape_calculator.CalculateMaterialSpectralShapes(parameters.cross_sections_.value());
    std::cout << "get material Spectral Shape" << std::endl;
    auto material_spectral_shapes = material_spectral_shape_calculator.material_spectral_shapes();

    for (auto& [id, vector] : material_spectral_shapes) {
      std::cout << "Material " << id << " spectral shape: ";
      for (auto& val : vector)
        std::cout << val << ",";
      std::cout << std::endl;
    }

    std::cout << "One group xsec" << std::endl;
    auto one_group_cross_sections = std::make_shared<data::cross_sections::CollapsedOneGroupCrossSections>(
        *parameters.cross_sections_.value(), material_spectral_shapes);
    acceleration::two_grid::spectral_shape::DomainSpectralShapes<dim> domain_spectral_shape_calculator;
    std::cout << "Domain spectral shape" << std::endl;
    auto group_to_domain_spectral_shape_map = domain_spectral_shape_calculator.CalculateDomainSpectralShapes(
        material_spectral_shapes, *domain_ptr);
    std::cout << "Flux corrector" << std::endl;
    auto flux_corrector = std::make_unique<acceleration::two_grid::FluxCorrector>(group_to_domain_spectral_shape_map);
    std::cout << "Isotropic Residual" << std::endl;
    auto domain_isotropic_residual_ptr = std::make_unique<calculator::residual::DomainIsotropicResidual<dim>>(
        std::make_unique<calculator::residual::CellIsotropicResidual<dim>>(parameters.cross_sections_.value(),
                                                                           finite_element_ptr), domain_ptr);
    auto two_grid_parameters{ parameters };
    two_grid_parameters.name = "Two-grid diffusion";
    two_grid_parameters.use_two_grid_ = false;
    two_grid_parameters.framework_level_ = 1;
    two_grid_parameters.output_filename_base = parameters.output_filename_base + "_two_grid";
    two_grid_parameters.cross_sections_ = one_group_cross_sections;
    two_grid_parameters.two_grid_data_.one_group_cross_sections_ptr_ = one_group_cross_sections;
    two_grid_parameters.equation_type = problem::EquationType::kTwoGridDiffusion;
    two_grid_parameters.eigen_solver_type = std::nullopt;
    two_grid_parameters.group_solver_type = problem::InGroupSolverType::kSourceIteration;
    two_grid_parameters.neutron_energy_groups = 1;

    std::unique_ptr<FrameworkI> subroutine_framework_ptr{ nullptr };
    if (subroutine_framework_helper_ptr_ != nullptr) {
      subroutine_framework_ptr = std::move(subroutine_framework_helper_ptr_->BuildFramework(builder, two_grid_parameters));
    } else {
      subroutine_framework_ptr = std::move(BuildFramework(builder, two_grid_parameters));
    }
    auto rhs_vector = std::make_shared<dealii::Vector<double>>(system_ptr->current_moments->GetMoment({0, 0, 0}).size());

    auto framework_ptr = dynamic_cast<Framework*>(subroutine_framework_ptr.get());
    auto outer_iteration_ptr = dynamic_cast<iteration::outer::OuterFixedSourceIteration*>(framework_ptr->outer_iterator_ptr());
    AssertThrow(outer_iteration_ptr != nullptr, dealii::ExcMessage("Error building Two-grid, outer iteration dynamic pointer null"))
    auto group_iteration_ptr = dynamic_cast<typename iteration::group::GroupSourceIteration<dim>*>(outer_iteration_ptr->group_iterator_ptr());
    AssertThrow(group_iteration_ptr != nullptr, dealii::ExcMessage("Error building Two-grid, group iteration dynamic pointer null"))
    auto fixed_upater_ptr = dynamic_cast<formulation::updater::FixedUpdater<dim>*>(group_iteration_ptr->source_updater_ptr());
    AssertThrow(fixed_upater_ptr != nullptr, dealii::ExcMessage("Error building Two-grid, fixed updater dynamic pointer null"))
    fixed_upater_ptr->SetRHSConstant(rhs_vector);

    group_post_processing_subroutine = std::make_unique<iteration::subroutine::TwoGridAcceleration>(
        std::move(flux_corrector), std::move(subroutine_framework_ptr),
        std::move(domain_isotropic_residual_ptr), rhs_vector);
    std::cout << "Two grid setup complete ==========================================================================\n";
  }

  if (group_post_processing_subroutine != nullptr)
    group_iteration_ptr->AddPostIterationSubroutine(std::move(group_post_processing_subroutine));

  std::unique_ptr<iteration::subroutine::SubroutineI> post_processing_subroutine{ nullptr };
  // For NDA we need to build a sub-routine framework to run the NDA process
  if (parameters.use_nda_) {
    auto nda_parameters{ parameters };
    nda_parameters.name = "NDA Drift-Diffusion";
    nda_parameters.use_nda_ = false;
    nda_parameters.framework_level_ = 1;
    nda_parameters.output_filename_base = parameters.output_filename_base + "_nda";
    nda_parameters.nda_data_.angular_flux_integrator_ptr_ =
        Shared(builder.BuildAngularFluxIntegrator(quadrature_set_ptr));
    nda_parameters.nda_data_.higher_order_moments_ptr_ = system_ptr->current_moments;
    nda_parameters.nda_data_.higher_order_angular_flux_ = angular_solutions_;
    std::unique_ptr<FrameworkI> subroutine_framework_ptr{ nullptr };
    if (subroutine_framework_helper_ptr_ != nullptr) {
      subroutine_framework_ptr = std::move(subroutine_framework_helper_ptr_->BuildFramework(builder, nda_parameters));
    } else {
      subroutine_framework_ptr = std::move(BuildFramework(builder, nda_parameters));
    }
    post_processing_subroutine = builder.BuildSubroutine(std::move(subroutine_framework_ptr),
                                                         iteration::subroutine::SubroutineName::kGetScalarFluxFromFramework);
  }

  std::unique_ptr<OuterIteration> outer_iteration_ptr{ nullptr };

  if (parameters.eigen_solver_type.has_value()){
    if (parameters.k_effective_updater == eigenvalue::k_eigenvalue::K_EffectiveUpdaterName::kCalculatorViaRayleighQuotient) {
      outer_iteration_ptr = builder.BuildOuterIteration(std::move(group_iteration_ptr),
                                                        builder.BuildParameterConvergenceChecker(1e-6, 1000),
                                                        builder.BuildKEffectiveUpdater(),
                                                        updater_pointers.fission_source_updater_ptr,
                                                        parameters.output_filename_base);
    } else {
      outer_iteration_ptr = builder.BuildOuterIteration(std::move(group_iteration_ptr),
                                                        builder.BuildParameterConvergenceChecker(1e-6, 1000),
                                                        builder.BuildKEffectiveUpdater(finite_element_ptr,
                                                                                       parameters.cross_sections_.value(),
                                                                                       domain_ptr),
                                                        updater_pointers.fission_source_updater_ptr,
                                                        parameters.output_filename_base);
    }
  } else {
    outer_iteration_ptr = builder.BuildOuterIteration(std::move(group_iteration_ptr),
                                                      builder.BuildParameterConvergenceChecker(1e-6, 1000),
                                                      parameters.output_filename_base);
  }

  // Add subroutines if applicable
  if (post_processing_subroutine != nullptr) {
    outer_iteration_ptr->AddPostIterationSubroutine(std::move(post_processing_subroutine));
  }

  if (parameters.output_scattering_source_as_vtu) {
    try {
      // Install if port is present
      auto vector_to_vtu_instrument = std::make_shared<instrumentation::BasicInstrument<dealii::Vector<double>>>(
          std::make_unique<typename instrumentation::outstream::VectorToVTU<dim>>(domain_ptr,
                                                                                  "scattering_source",
                                                                                  "scattering_source",
                                                                                  parameters.output_filename_base
                                                                                      + "_scattering_source"));
      instrumentation::GetPort<iteration::outer::data_names::ScatteringSourcePort>(*outer_iteration_ptr)
          .AddInstrument(vector_to_vtu_instrument);
    } catch (std::bad_cast &) {
      AssertThrow(false, dealii::ExcMessage("Error installing scattering source to vtu instrument, port is not present"))
    }
  }
  if (parameters.output_fission_source_as_vtu) {
    try {
      auto fission_source_vector_to_vtu_instrument =
          std::make_shared<instrumentation::BasicInstrument<dealii::Vector<double>>>(
              std::make_unique<typename instrumentation::outstream::VectorToVTU<dim>>(domain_ptr,
                                                                                      "fission_source",
                                                                                      "fission_source",
                                                                                      parameters.output_filename_base
                                                                                          + "_fission_source"));
      instrumentation::GetPort<iteration::outer::data_names::FissionSourcePort>(*outer_iteration_ptr)
          .AddInstrument(fission_source_vector_to_vtu_instrument);
    } catch (std::bad_cast &) {
      AssertThrow(false, dealii::ExcMessage("Error installing fission source  to vtu instrument, port is not present"))
    }
  }
  if (parameters.output_scalar_flux_as_vtu) {
    try {
      auto scalar_flux_to_vtu_instrument = std::make_shared<instrumentation::BasicInstrument<dealii::Vector<double>>>(
          std::make_unique<typename instrumentation::outstream::VectorToVTU<dim>>(domain_ptr,
                                                                                  "scalar_flux",
                                                                                  "scalar_flux",
                                                                                  parameters.output_filename_base
                                                                                      + "_scalar_flux"));
      instrumentation::GetPort<iteration::outer::data_names::ScalarFluxPort>(*outer_iteration_ptr)
          .AddInstrument(scalar_flux_to_vtu_instrument);
    } catch (std::bad_cast &) {
      AssertThrow(false, dealii::ExcMessage("Error scalar flux to vtu instrument, port is not present"))
    }
  }

  validator.ReportValidation();

  return std::make_unique<framework::Framework>(std::move(system_ptr),
                                                std::move(initializer_ptr),
                                                std::move(outer_iteration_ptr),
                                                std::make_unique<results::OutputDealiiVtu<dim>>(domain_ptr));

}

template<int dim>
auto FrameworkHelper<dim>::BuildFramework(
    builder::FrameworkBuilderI<dim>& builder,
    FrameworkParameters& parameters,
    system::moments::SphericalHarmonicI* previous_solution_ptr) -> std::unique_ptr<framework::FrameworkI> {
  auto framework_ptr = BuildFramework(builder, parameters);
  auto dynamic_framework_ptr = dynamic_cast<framework::Framework*>(framework_ptr.get());
  auto outer_iteration_ptr = dynamic_framework_ptr->outer_iterator_ptr();

  auto fourier_instrument = Shared(
      InstrumentBuilder::BuildInstrument<system::moments::SphericalHarmonicI>(
          InstrumentName::kFourierTransformOfAllGroupScalarFluxErrorToFile,
          previous_solution_ptr,
          parameters.output_filename_base + "_fourier_of_error_group"));
  using FourierDataPort = iteration::outer::data_names::SolutionMomentsPort;
  instrumentation::GetPort<FourierDataPort>(*outer_iteration_ptr).AddInstrument(fourier_instrument);
  return framework_ptr;
}

template class FrameworkHelper<1>;
template class FrameworkHelper<2>;
template class FrameworkHelper<3>;

} // namespace bart::framework
