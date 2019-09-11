#include "framework/builder/cfem_framework_builder.h"

#include <fstream>
#include <streambuf>
#include <deal.II/base/mpi.h>

#include "calculator/cell/integrated_fission_source.h"
#include "calculator/cell/total_aggregated_fission_source.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/reporter/mpi_noisy.h"
#include "domain/definition.h"
#include "domain/finite_element/finite_element_gaussian.h"
#include "domain/mesh/mesh_cartesian.h"
#include "eigenvalue/k_effective/updater_via_fission_source.h"
#include "formulation/cfem_diffusion_stamper.h"
#include "formulation/scalar/cfem_diffusion.h"
#include "framework/framework.h"
#include "iteration/updater/source_updater_gauss_seidel.h"
#include "iteration/updater/fixed_updater.h"
#include "iteration/initializer/set_fixed_terms_once.h"
#include "iteration/group/group_source_iteration.h"
#include "iteration/outer/outer_power_iteration.h"
#include "material/material_protobuf.h"
#include "problem/parameter_types.h"
#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/final_checker_or_n.h"
#include "quadrature/angular/angular_quadrature_scalar.h"
#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "solver/group/single_group_solver.h"
#include "solver/gmres.h"
#include "system/system.h"
#include "system/solution/mpi_group_angular_solution.h"
#include "system/terms/term.h"
#include "system/terms/term_types.h"
#include "system/moments/spherical_harmonic.h"


namespace bart {

namespace framework {

namespace builder {

template <int dim>
std::unique_ptr<FrameworkI> CFEM_FrameworkBuilder<dim>::BuildFramework(
    problem::ParametersI &prm,
    dealii::ParameterHandler &d2_prm) {

  std::cout << "Setting up materials" << std::endl;

  std::ifstream mapping_file(prm.MaterialMapFilename());
  std::string material_mapping(
      (std::istreambuf_iterator<char>(mapping_file)),
      std::istreambuf_iterator<char>());

  const int n_groups = 1;

  MaterialProtobuf materials(d2_prm);
  auto cross_sections_ptr = std::make_shared<bart::data::CrossSections>(materials);

  std::cout << "Building Finite Element"<< std::endl;
  std::shared_ptr<FiniteElement>finite_element_ptr(std::move(BuildFiniteElement(
      &prm)));

  std::cout << "Building Domain" << std::endl;
  std::shared_ptr<Domain> domain_ptr(std::move(BuildDomain(
      &prm,finite_element_ptr, material_mapping)));

  std::cout << "Setting up domain" << std::endl;
  domain_ptr->SetUpMesh().SetUpDOF();

  std::cout << "Building Stamper" << std::endl;
  std::shared_ptr<CFEMStamper> stamper_ptr(std::move(BuildStamper(
      &prm, domain_ptr, finite_element_ptr, cross_sections_ptr)));

  std::cout << "Building fixed updater" << std::endl;
  auto fixed_updater_ptr = BuildFixedUpdater(stamper_ptr);

  std::cout << "Building source updater" << std::endl;
  std::shared_ptr<SourceUpdater> source_updater_ptr(
      std::move(BuildSourceUpdater(&prm, stamper_ptr)));

  std::cout << "Building Initializer" << std::endl;
  auto initializer_ptr =
      std::make_unique<iteration::initializer::SetFixedTermsOnce>(
          std::move(fixed_updater_ptr), prm.NEnergyGroups(), n_groups);

  std::cout << "Building single group solver" << std::endl;
  auto single_group_solver_ptr = BuildSingleGroupSolver();


  std::cout << "Building inner iteration objects" << std::endl;

  auto in_group_final_checker = BuildMomentConvergenceChecker(1e-10, 100);

  // Build reporter
  using Reporter = bart::convergence::reporter::MpiNoisy;
  int this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  auto pout_ptr = std::make_unique<dealii::ConditionalOStream>(std::cout, this_process == 0);
  auto reporter = std::make_shared<Reporter>(std::move(pout_ptr));

  // Scalar Quadrature
  using ScalarQuadrature = quadrature::angular::AngularQuadratureScalar<dim>;
  auto quadrature_ptr = std::make_shared<ScalarQuadrature>();

  // Moment calculator
  using MomentCalculator = quadrature::calculators::SphericalHarmonicZerothMoment<dim>;
  auto moment_calculator = std::make_unique<MomentCalculator>(quadrature_ptr);

  // Solution group
  auto solution_ptr =
      std::make_shared<system::solution::MPIGroupAngularSolution>(n_groups);

  data::MatrixParameters matrix_param;
  domain_ptr->FillMatrixParameters(matrix_param, prm.Discretization());

  std::cout << "==Filling solution object" << std::endl;

  auto& solution = solution_ptr->operator[](0);
  solution.reinit(matrix_param.rows, MPI_COMM_WORLD);
  auto local_elements = solution.locally_owned_elements();
  for (auto index : local_elements) {
    solution[index] = 1.0;
  }
  solution.compress(dealii::VectorOperation::insert);

  std::cout << "Building in-group iteration" << std::endl;

  auto in_group_iteration =
      std::make_unique<iteration::group::GroupSourceIteration<dim>>(
          std::move(single_group_solver_ptr),
          std::move(in_group_final_checker),
          std::move(moment_calculator),
          solution_ptr,
          source_updater_ptr,
          reporter);

  std::cout << "Building K_Effective updater" << std::endl;

  // KEffectiveUpdater
  using FissionSourceCalulator = calculator::cell::TotalAggregatedFissionSource<dim>;
  using IntegratedFissionSourceCalc = calculator::cell::IntegratedFissionSource<dim>;
  auto int_fission_ptr = std::make_unique<IntegratedFissionSourceCalc>(
      finite_element_ptr, cross_sections_ptr);
  auto fission_source_calc_ptr = std::make_unique<FissionSourceCalulator>(
      std::move(int_fission_ptr), domain_ptr);
  using KEffectiveUpdater = eigenvalue::k_effective::UpdaterViaFissionSource;

  // KEffective convergence checker
  using KEffConvChecker = bart::convergence::parameters::SingleParameterChecker;
  auto k_eff_conv_checker = std::make_unique<KEffConvChecker>();
  using KEffectiveFinalCovergence = bart::convergence::FinalCheckerOrN<double,
      bart::convergence::parameters::SingleParameterChecker>;
  auto k_eff_final_checker = std::make_unique<KEffectiveFinalCovergence>(
      std::move(k_eff_conv_checker)
      );

  // KEffective updater
  auto k_effective_updater = std::make_unique<KEffectiveUpdater>(
      std::move(fission_source_calc_ptr), 2.0, 10);

  std::cout << "Building outer-iteration" << std::endl;

  using PowerIteration = iteration::outer::OuterPowerIteration;
  auto power_iteration_ptr = std::make_unique<PowerIteration>(
      std::move(in_group_iteration),
      std::move(k_eff_final_checker),
      std::move(k_effective_updater),
      source_updater_ptr,
      reporter
  );

  std::cout << "Building System" << std::endl;

  auto system = std::make_unique<system::System>();

  system->total_groups = prm.NEnergyGroups();
  std::unordered_set<bart::system::terms::VariableLinearTerms>
      source_terms{bart::system::terms::VariableLinearTerms::kScatteringSource,
                   bart::system::terms::VariableLinearTerms::kFissionSource};
  system->right_hand_side_ptr_ =
      std::make_unique<system::terms::MPILinearTerm>(source_terms);
  system->left_hand_side_ptr_ =
      std::make_unique<system::terms::MPIBilinearTerm>();

  std::cout << "Filling system" << std::endl;

  // Fill system with objects
  for (int group = 0; group <= prm.NEnergyGroups(); ++group) {

    // LHS
    auto fixed_matrix_ptr = bart::data::BuildMatrix(matrix_param);
    system->left_hand_side_ptr_->SetFixedTermPtr(group, fixed_matrix_ptr);

    // RHS
    auto fixed_vector_ptr =
        std::make_shared<bart::system::MPIVector>(matrix_param.rows, MPI_COMM_WORLD);
    system->right_hand_side_ptr_->SetFixedTermPtr(group, fixed_vector_ptr);

    for (auto term : source_terms) {
      auto variable_vector_ptr =
          std::make_shared<bart::system::MPIVector>(matrix_param.rows, MPI_COMM_WORLD);
      system->right_hand_side_ptr_->SetVariableTermPtr(
          group, term, variable_vector_ptr);
    }
    // Moments
    system->current_moments =
        std::make_unique<system::moments::SphericalHarmonic>(prm.NEnergyGroups(), 0);
    system->previous_moments =
        std::make_unique<system::moments::SphericalHarmonic>(prm.NEnergyGroups(), 0);
  }

  std::cout << "Fill system moments" << std::endl;

  for (auto& moment_pair : system->current_moments->moments()) {
    auto index = moment_pair.first;
    auto& current_moment = system->current_moments->operator[](index);
    current_moment.reinit(solution_ptr->operator[](0).size());
    auto& previous_moment = system->previous_moments->operator[](index);
    previous_moment.reinit(solution_ptr->operator[](0).size());
    current_moment = 1;
    previous_moment = 1;
  }

  // Initialize System
  system->k_effective = 1.16;
  system->total_groups = prm.NEnergyGroups();
  system->total_angles = 1;

  return std::make_unique<framework::Framework>(
      std::move(system),
      std::move(initializer_ptr),
      std::move(power_iteration_ptr));
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildFiniteElement(
    problem::ParametersI *problem_parameters)-> std::unique_ptr<FiniteElement> {
  return std::make_unique<domain::finite_element::FiniteElementGaussian<dim>>(
      problem::DiscretizationType::kContinuousFEM,
      problem_parameters->FEPolynomialDegree());
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildDomain(
    problem::ParametersI *problem_parameters,
    const std::shared_ptr<FiniteElement> &finite_element_ptr,
    std::string material_mapping)-> std::unique_ptr<Domain> {

  // Build mesh
  auto mesh_ptr = std::make_unique<domain::mesh::MeshCartesian<dim>>(
      problem_parameters->SpatialMax(),
          problem_parameters->NCells(),
          material_mapping);

  return std::make_unique<domain::Definition<dim>>(
      std::move(mesh_ptr), finite_element_ptr);
}
template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildStamper(
    problem::ParametersI *problem_parameters,
    const std::shared_ptr<Domain> &domain_ptr,
    const std::shared_ptr<FiniteElement> &finite_element_ptr,
    const std::shared_ptr<CrossSections> &cross_sections_ptr)
-> std::unique_ptr<CFEMStamper> {

  std::unique_ptr<CFEMStamper> return_ptr = nullptr;

  // Diffusion Stamper
  if (problem_parameters->TransportModel() == problem::EquationType::kDiffusion) {

    auto diffusion_ptr = std::make_unique<formulation::scalar::CFEM_Diffusion<dim>>(
        finite_element_ptr, cross_sections_ptr);

    return_ptr = std::move(
        std::make_unique<formulation::CFEM_DiffusionStamper<dim>>(
            std::move(diffusion_ptr),
            domain_ptr,
            problem_parameters->ReflectiveBoundary()));

  } else {
    AssertThrow(false, dealii::ExcMessage("Unsuppored equation type passed"
                                         "to BuildScalarFormulation"));
  }

  return return_ptr;
}
template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildSourceUpdater(
    problem::ParametersI *,
    const std::shared_ptr<CFEMStamper> stamper_ptr)
    -> std::unique_ptr<SourceUpdater> {
  // TODO(Josh): Add option for non-gauss-seidel updating
  using SourceUpdater = iteration::updater::SourceUpdaterGaussSeidel<CFEMStamper>;
  return std::make_unique<SourceUpdater>(stamper_ptr);
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildParameterConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<ParameterConvergenceChecker>{

  using CheckerType = convergence::parameters::SingleParameterChecker;
  using FinalCheckerType = convergence::FinalCheckerOrN<double, CheckerType>;

  auto single_checker_ptr = std::make_unique<CheckerType>(max_delta);
  auto return_ptr = std::make_unique<FinalCheckerType>(
      std::move(single_checker_ptr));

  return_ptr->SetMaxIterations(max_iterations);

  return std::move(return_ptr);
}

template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildMomentConvergenceChecker(
    double max_delta, int max_iterations)
-> std::unique_ptr<MomentConvergenceChecker>{
  //TODO(Josh): Add option for using other than L1Norm

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
auto CFEM_FrameworkBuilder<dim>::BuildFixedUpdater(
    const std::shared_ptr<CFEMStamper> &stamper_ptr)
-> std::unique_ptr<CFEM_FrameworkBuilder::FixedUpdater> {
  using FixedUpdaterType = iteration::updater::FixedUpdater<CFEMStamper>;
  return std::make_unique<FixedUpdaterType>(stamper_ptr);
}
template<int dim>
auto CFEM_FrameworkBuilder<dim>::BuildSingleGroupSolver(
    const int max_iterations,
    const double convergence_tolerance) -> std::unique_ptr<SingleGroupSolver> {
  std::unique_ptr<SingleGroupSolver> return_ptr = nullptr;

  auto linear_solver_ptr = std::make_unique<solver::GMRES>(max_iterations,
                                                    convergence_tolerance);

  return_ptr = std::move(
      std::make_unique<solver::group::SingleGroupSolver>(
          std::move(linear_solver_ptr)));

  return return_ptr;
}

template class CFEM_FrameworkBuilder<1>;
template class CFEM_FrameworkBuilder<2>;
template class CFEM_FrameworkBuilder<3>;



} // namespace builder

} // namespace framework

} // namespace bart