#include <eigenvalue/k_effective/updater_via_rayleigh_quotient.hpp>
#include "framework/framework_helper.hpp"

#include "convergence/tests/final_checker_mock.h"
#include "eigenvalue/k_effective/tests/k_effective_updater_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/updater/tests/boundary_conditions_updater_mock.h"
#include "formulation/updater/tests/fission_source_updater_mock.h"
#include "formulation/updater/tests/fixed_updater_mock.h"
#include "formulation/updater/tests/scattering_source_updater_mock.h"
#include "formulation/scalar/tests/diffusion_mock.h"
#include "framework/builder/framework_builder_i.hpp"
#include "framework/builder/tests/framework_builder_mock.hpp"
#include "framework/builder/tests/framework_validator_mock.hpp"
#include "framework/framework_parameters.hpp"
#include "framework/framework.hpp"
#include "iteration/initializer/tests/initializer_mock.h"
#include "iteration/group/tests/group_solve_iteration_mock.h"
#include "iteration/outer/tests/outer_iteration_mock.hpp"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "domain/tests/definition_mock.h"
#include "material/tests/mock_material.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/calculators/tests/spherical_harmonic_moments_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "solver/group/tests/single_group_solver_mock.h"
#include "system/tests/system_helper_mock.hpp"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/system.hpp"

namespace  {

using namespace bart;
using ::testing::Return, ::testing::ByMove, ::testing::DoDefault, ::testing::_, ::testing::NiceMock;
using ::testing::Ref, ::testing::Pointee, ::testing::ReturnRef, ::testing::ContainerEq, ::testing::SizeIs;
using ::testing::A, ::testing::AllOf;
using ::testing::NotNull, ::testing::WhenDynamicCastTo;
using ::testing::AtLeast;

using FrameworkPart = framework::builder::FrameworkPart;

template <typename DimensionWrapper>
class FrameworkHelperBuildFrameworkIntegrationTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  // Mock types
  using DiffusionFormulationMock = typename formulation::scalar::DiffusionMock<dim>;
  using DomainMock = typename domain::DefinitionMock<dim>;
  using FiniteElementMock = typename domain::finite_element::FiniteElementMock<dim>;
  using FrameworkBuidler = framework::builder::FrameworkBuilderMock<dim>;
  using FrameworkParameters = framework::FrameworkParameters;
  using GroupSolutionMock = system::solution::MPIGroupAngularSolutionMock;
  using GroupSolveIterationMock = iteration::group::GroupSolveIterationMock;
  using InitializerMock = iteration::initializer::InitializerMock;
  using KEffectiveUpdaterMock = eigenvalue::k_effective::K_EffectiveUpdaterMock;
  using MomentCalculatorMock = quadrature::calculators::SphericalHarmonicMomentsMock;
  using MomentConvergenceCheckerMock = convergence::FinalCheckerMock<system::moments::MomentVector>;
  using MomentMapConvergenceCheckerMock = convergence::FinalCheckerMock<const system::moments::MomentsMap>;
  using OuterIterationMock = iteration::outer::OuterIterationMock;
  using ParameterConvergenceCheckerMock = convergence::FinalCheckerMock<double>;
  using QuadratureSetMock = typename quadrature::QuadratureSetMock<dim>;
  using StamperMock = typename formulation::StamperMock<dim>;
  using SAAFFormulationMock = typename formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using SingleGroupSolverMock = solver::group::SingleGroupSolverMock;
  using SystemHelper = typename system::SystemHelperMock<dim>;
  using System = system::System;
  using Validator = framework::builder::FrameworkValidatorMock;

  using UpdaterPointers = FrameworkBuidler::UpdaterPointers;
  using BoundaryConditionsUpdaterMock = formulation::updater::BoundaryConditionsUpdaterMock;
  using FissionSourceUpdaterMock = formulation::updater::FissionSourceUpdaterMock;
  using FixedTermUpdaterMock = formulation::updater::FixedUpdaterMock;
  using ScatteringSourceUpdaterMock = formulation::updater::ScatteringSourceUpdaterMock;

  // Test Type and object
  using FrameworkHelper = typename framework::FrameworkHelper<dim>;
  std::unique_ptr<FrameworkHelper> test_helper_ptr_;

  // Mock pointers and observation pointers
  DiffusionFormulationMock* diffusion_formulation_obs_ptr_{ nullptr };
  DomainMock* domain_obs_ptr_{ nullptr };
  FiniteElementMock* finite_element_obs_ptr_{ nullptr };
  GroupSolutionMock* group_solution_obs_ptr_{ nullptr };
  GroupSolveIterationMock* group_solve_iteration_obs_ptr{ nullptr };
  InitializerMock* initializer_obs_ptr_{ nullptr };
  KEffectiveUpdaterMock* k_effective_updater_obs_ptr_{ nullptr };
  KEffectiveUpdaterMock* k_effective_updater_rayleigh_obs_ptr_{ nullptr };
  MomentCalculatorMock* moment_calculator_obs_ptr_{ nullptr };
  MomentConvergenceCheckerMock* moment_convergence_checker_obs_ptr_{ nullptr };
  MomentMapConvergenceCheckerMock* moment_map_convergence_checker_obs_ptr_{ nullptr };
  OuterIterationMock* outer_iteration_obs_ptr_{ nullptr };
  OuterIterationMock* outer_iteration_eigensolve_obs_ptr_{ nullptr };
  ParameterConvergenceCheckerMock* parameter_convergence_checker_obs_ptr_{ nullptr };
  std::shared_ptr<QuadratureSetMock> quadrature_set_mock_ptr_{ nullptr };
  SAAFFormulationMock* saaf_formulation_obs_ptr_{ nullptr };
  SingleGroupSolverMock* single_group_solver_obs_ptr_{ nullptr };
  StamperMock* stamper_obs_ptr_{ nullptr };
  std::shared_ptr<SystemHelper> system_helper_mock_ptr_{ nullptr };
  System* system_obs_ptr_{ nullptr };

  UpdaterPointers updater_pointers_;

  // Test parameters and supporting objects
  FrameworkBuidler mock_builder_;
  FrameworkParameters default_parameters_;
  const int total_quadrature_angles{ test_helpers::RandomInt(10, 20) };
  std::vector<domain::CellPtr<dim>> cells_;
  system::MPIVector solution_;
  Validator mock_validator_;

  auto RunTest(const FrameworkParameters& parameters) -> void;
  auto SetUp() -> void override;
};

template <typename T>
auto ReturnByMove(T& to_return) {
  return Return(ByMove(std::move(to_return)));
}

template <typename DimensionWrapper>
auto FrameworkHelperBuildFrameworkIntegrationTests<DimensionWrapper>::SetUp() -> void {
  using DomainSize = FrameworkParameters::DomainSize;
  using NumberOfCells = FrameworkParameters::NumberOfCells;
  NiceMock<btest::MockMaterial> mock_material;

  // Default framework parameters (for tests)
  default_parameters_.neutron_energy_groups = test_helpers::RandomInt(1, 4);
  default_parameters_.cross_sections_ = std::make_shared<data::CrossSections>(mock_material);
  default_parameters_.material_mapping = "1 1";
  default_parameters_.domain_size = DomainSize(test_helpers::RandomVector(dim, 0, 100));
  default_parameters_.uniform_refinements = test_helpers::RandomInt(1, 4);

  auto random_n_cells_double = test_helpers::RandomVector(dim, 10, 20);
  std::vector<int> random_n_cells(random_n_cells_double.cbegin(), random_n_cells_double.cend());
  default_parameters_.number_of_cells = NumberOfCells(random_n_cells);
  int total_cells = std::accumulate(random_n_cells.cbegin(), random_n_cells.cend(), 1.0,
                                    std::multiplies<int>());
  for (int i = 0; i < total_cells; ++i) {
    cells_.emplace_back();
  }

  // Mocks and observation pointers
  auto diffusion_formulation_ptr = std::make_unique<NiceMock<DiffusionFormulationMock>>();
  diffusion_formulation_obs_ptr_ = diffusion_formulation_ptr.get();
  auto domain_ptr = std::make_unique<NiceMock<DomainMock>>();
  domain_obs_ptr_ = domain_ptr.get();
  auto finite_element_ptr = std::make_unique<NiceMock<FiniteElementMock>>();
  finite_element_obs_ptr_ = finite_element_ptr.get();
  auto group_solution_ptr = std::make_unique<NiceMock<GroupSolutionMock>>();
  group_solution_obs_ptr_ = group_solution_ptr.get();
  auto group_solve_iteration_ptr = std::make_unique<NiceMock<GroupSolveIterationMock>>();
  group_solve_iteration_obs_ptr = group_solve_iteration_ptr.get();
  auto initializer_ptr = std::make_unique<NiceMock<InitializerMock>>();
  initializer_obs_ptr_ = initializer_ptr.get();
  auto k_effective_updater_ptr = std::make_unique<NiceMock<KEffectiveUpdaterMock>>();
  k_effective_updater_obs_ptr_ = k_effective_updater_ptr.get();
  auto k_effective_updater_rayleigh_ptr = std::make_unique<NiceMock<KEffectiveUpdaterMock>>();
  k_effective_updater_rayleigh_obs_ptr_ = k_effective_updater_rayleigh_ptr.get();
  auto moment_calculator_ptr = std::make_unique<NiceMock<MomentCalculatorMock>>();
  moment_calculator_obs_ptr_ = moment_calculator_ptr.get();
  auto moment_convergence_checker_ptr = std::make_unique<NiceMock<MomentConvergenceCheckerMock>>();
  moment_convergence_checker_obs_ptr_ = moment_convergence_checker_ptr.get();
  auto moment_map_convergence_checker_ptr = std::make_unique<NiceMock<MomentMapConvergenceCheckerMock>>();
  moment_map_convergence_checker_obs_ptr_ = moment_map_convergence_checker_ptr.get();
  auto outer_iteration_ptr = std::make_unique<NiceMock<OuterIterationMock>>();
  outer_iteration_obs_ptr_ = outer_iteration_ptr.get();
  auto outer_iteration_eigensolve_ptr = std::make_unique<NiceMock<OuterIterationMock>>();
  outer_iteration_eigensolve_obs_ptr_ = outer_iteration_eigensolve_ptr.get();
  auto parameter_convergence_checker_ptr = std::make_unique<NiceMock<ParameterConvergenceCheckerMock>>();
  parameter_convergence_checker_obs_ptr_ = parameter_convergence_checker_ptr.get();
  quadrature_set_mock_ptr_ = std::make_shared<QuadratureSetMock>();
  auto saaf_ptr = std::make_unique<NiceMock<SAAFFormulationMock>>();
  saaf_formulation_obs_ptr_ = saaf_ptr.get();
  auto single_group_solver_ptr = std::make_unique<NiceMock<SingleGroupSolverMock>>();
  single_group_solver_obs_ptr_ = single_group_solver_ptr.get();
  auto stamper_ptr = std::make_unique<NiceMock<StamperMock>>();
  stamper_obs_ptr_ = stamper_ptr.get();
  system_helper_mock_ptr_ = std::make_shared<NiceMock<SystemHelper>>();
  auto system_ptr = std::make_unique<System>();
  system_obs_ptr_ = system_ptr.get();

  updater_pointers_.boundary_conditions_updater_ptr = std::make_shared<NiceMock<BoundaryConditionsUpdaterMock>>();
  updater_pointers_.fission_source_updater_ptr = std::make_shared<NiceMock<FissionSourceUpdaterMock>>();
  updater_pointers_.fixed_updater_ptr = std::make_shared<NiceMock<FixedTermUpdaterMock>>();
  updater_pointers_.scattering_source_updater_ptr = std::make_shared<NiceMock<ScatteringSourceUpdaterMock>>();


  ON_CALL(*group_solution_obs_ptr_, GetSolution(_)).WillByDefault(ReturnRef(solution_));

  using DiffusionFormulationPtr = std::unique_ptr<typename FrameworkBuidler::DiffusionFormulation>;
  using SAAFFormulationPtr = std::unique_ptr<typename FrameworkBuidler::SAAFFormulation>;

  ON_CALL(mock_builder_, BuildDiffusionFormulation(_,_,_)).WillByDefault(ReturnByMove(diffusion_formulation_ptr));
  ON_CALL(mock_builder_, BuildDomain(_, _, _, _)).WillByDefault(ReturnByMove(domain_ptr));
  ON_CALL(mock_builder_, BuildFiniteElement(_,_,_)).WillByDefault(ReturnByMove(finite_element_ptr));
  ON_CALL(mock_builder_, BuildGroupSolution(_)).WillByDefault(ReturnByMove(group_solution_ptr));
  ON_CALL(mock_builder_, BuildGroupSolveIteration(_,_,_,_,_,_)).WillByDefault(ReturnByMove(group_solve_iteration_ptr));
  ON_CALL(mock_builder_, BuildInitializer(_,_,_)).WillByDefault(ReturnByMove(initializer_ptr));
  ON_CALL(mock_builder_, BuildKEffectiveUpdater(_,_,_)).WillByDefault(ReturnByMove(k_effective_updater_ptr));
  ON_CALL(mock_builder_, BuildKEffectiveUpdater()).WillByDefault(ReturnByMove(k_effective_updater_rayleigh_ptr));
  ON_CALL(mock_builder_, BuildMomentCalculator(_)).WillByDefault(ReturnByMove(moment_calculator_ptr));
  ON_CALL(mock_builder_, BuildMomentCalculator(_,_)).WillByDefault(ReturnByMove(moment_calculator_ptr));
  ON_CALL(mock_builder_, BuildMomentConvergenceChecker(_,_)).WillByDefault(ReturnByMove(moment_convergence_checker_ptr));
  ON_CALL(mock_builder_, BuildMomentMapConvergenceChecker(_,_))
      .WillByDefault(ReturnByMove(moment_map_convergence_checker_ptr));
  ON_CALL(mock_builder_, BuildOuterIteration(_,_,_)).WillByDefault(ReturnByMove(outer_iteration_ptr));
  ON_CALL(mock_builder_, BuildOuterIteration(_,_,_,_,_)).WillByDefault(ReturnByMove(outer_iteration_eigensolve_ptr));
  ON_CALL(mock_builder_, BuildParameterConvergenceChecker(_,_))
      .WillByDefault(ReturnByMove(parameter_convergence_checker_ptr));
  ON_CALL(mock_builder_, BuildQuadratureSet(_,_)).WillByDefault(Return(quadrature_set_mock_ptr_));
  ON_CALL(mock_builder_, BuildSAAFFormulation(_,_,_,_)).WillByDefault(ReturnByMove(saaf_ptr));
  ON_CALL(mock_builder_, BuildSingleGroupSolver(_,_)).WillByDefault(ReturnByMove(single_group_solver_ptr));
  ON_CALL(mock_builder_, BuildStamper(_)).WillByDefault(ReturnByMove(stamper_ptr));
  ON_CALL(mock_builder_, BuildUpdaterPointers(A<SAAFFormulationPtr>(),_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildUpdaterPointers(A<DiffusionFormulationPtr>(),_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildUpdaterPointers(_,_,_,_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildSystem(_,_,_,_,_,_)).WillByDefault(ReturnByMove(system_ptr));
  ON_CALL(mock_builder_, set_color_status_instrument_ptr(_)).WillByDefault(ReturnRef(mock_builder_));
  ON_CALL(mock_builder_, set_convergence_status_instrument_ptr(_)).WillByDefault(ReturnRef(mock_builder_));
  ON_CALL(mock_builder_, set_status_instrument_ptr(_)).WillByDefault(ReturnRef(mock_builder_));
  ON_CALL(mock_builder_, validator_ptr()).WillByDefault(Return(&mock_validator_));

  ON_CALL(*group_solve_iteration_obs_ptr, UpdateThisAngularSolutionMap(_))
      .WillByDefault(ReturnRef(*group_solve_iteration_obs_ptr));

  ON_CALL(*domain_obs_ptr_, SetUpMesh(_)).WillByDefault(ReturnRef(*domain_obs_ptr_));
  ON_CALL(*domain_obs_ptr_, SetUpDOF()).WillByDefault(ReturnRef(*domain_obs_ptr_));
  ON_CALL(*domain_obs_ptr_, Cells()).WillByDefault(Return(cells_));

  ON_CALL(*quadrature_set_mock_ptr_, size()).WillByDefault(Return(this->total_quadrature_angles));

  test_helper_ptr_ = std::make_unique<FrameworkHelper>(system_helper_mock_ptr_);
}

template<typename DimensionWrapper>
auto FrameworkHelperBuildFrameworkIntegrationTests<DimensionWrapper>::RunTest(
    const FrameworkParameters &parameters) -> void {
  auto& mock_builder = this->mock_builder_;
  int n_angles{ 1 };
  bool need_angular_storage{ false };
  const bool is_eigenvalue_solve {parameters.eigen_solver_type.has_value() };

  EXPECT_CALL(mock_builder, validator_ptr()).Times(AtLeast(1)).WillRepeatedly(DoDefault());
  EXPECT_CALL(mock_validator_, Parse(A<FrameworkParameters>()));


  // Mock Builder calls
  EXPECT_CALL(mock_builder, set_color_status_instrument_ptr(NotNull())).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, set_convergence_status_instrument_ptr(NotNull())).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, set_status_instrument_ptr(NotNull())).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildFiniteElement(parameters.cell_finite_element_type,
                                               parameters.discretization_type,
                                               parameters.polynomial_degree))
      .WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildDomain(parameters.domain_size,
                                        parameters.number_of_cells,
                                        Pointee(Ref(*this->finite_element_obs_ptr_)),
                                        parameters.material_mapping))
      .WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildStamper(Pointee(Ref(*this->domain_obs_ptr_))))
      .WillOnce(DoDefault());

  // Domain calls
  EXPECT_CALL(*this->domain_obs_ptr_, SetUpMesh(parameters.uniform_refinements)).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_obs_ptr_, SetUpDOF()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_obs_ptr_, Cells()).WillOnce(DoDefault());

  // Angular calls
  if (parameters.equation_type == problem::EquationType::kSelfAdjointAngularFlux) {
    // Angular equation types
    EXPECT_CALL(mock_builder, BuildQuadratureSet(parameters.angular_quadrature_type,
                                                 parameters.angular_quadrature_order.value()))
        .WillOnce(DoDefault());
    EXPECT_CALL(*quadrature_set_mock_ptr_, size()).WillOnce(DoDefault());
    EXPECT_CALL(mock_builder, BuildMomentCalculator(Pointee(Ref(*quadrature_set_mock_ptr_)),_)).WillOnce(DoDefault());
    n_angles = this->total_quadrature_angles;
  } else {
    // Scalar types
    EXPECT_CALL(mock_builder, BuildMomentCalculator(_)).WillOnce(DoDefault());
  }

  // Build a mapping of reflective boundaries
  std::map<problem::Boundary, bool> reflective_boundaries {
      {problem::Boundary::kXMin, false}, {problem::Boundary::kXMax, false},
      {problem::Boundary::kYMin, false}, {problem::Boundary::kYMax, false},
      {problem::Boundary::kZMin, false}, {problem::Boundary::kZMax, false},
  };
  for (auto& boundary : parameters.reflective_boundaries) {
    reflective_boundaries.at(boundary) = true;
  }

  // SAAF Specific calls
  if (parameters.equation_type == problem::EquationType::kSelfAdjointAngularFlux) {

    EXPECT_CALL(mock_builder, BuildSAAFFormulation(Pointee(Ref(*this->finite_element_obs_ptr_)),
                                                   Pointee(Ref(*parameters.cross_sections_.value())),
                                                   Pointee(Ref(*quadrature_set_mock_ptr_)),
                                                   formulation::SAAFFormulationImpl::kDefault))
        .WillOnce(DoDefault());
    EXPECT_CALL(*saaf_formulation_obs_ptr_, Initialize(cells_.at(0)));
    if (parameters.reflective_boundaries.empty()) {
      using SAAFFormulationPtr = std::unique_ptr<typename framework::builder::FrameworkBuilderI<dim>::SAAFFormulation>;
      EXPECT_CALL(mock_builder, BuildUpdaterPointers(
          A<SAAFFormulationPtr>(),
          Pointee(Ref(*stamper_obs_ptr_)),
          Pointee(Ref(*quadrature_set_mock_ptr_))))
          .WillOnce(DoDefault());
    } else {
      // Should have angular storage
      need_angular_storage = true;
      EXPECT_CALL(*system_helper_mock_ptr_, SetUpEnergyGroupToAngularSolutionPtrMap(
          _, parameters.neutron_energy_groups, total_quadrature_angles));
      EXPECT_CALL(mock_builder, BuildUpdaterPointers(
          Pointee(Ref(*saaf_formulation_obs_ptr_)),
          Pointee(Ref(*stamper_obs_ptr_)),
          Pointee(Ref(*quadrature_set_mock_ptr_)),
          ContainerEq(reflective_boundaries),
          _))
          .WillOnce(DoDefault());
      EXPECT_CALL(*group_solve_iteration_obs_ptr, UpdateThisAngularSolutionMap(_)).WillOnce(DoDefault());
    }
  } else if (parameters.equation_type == problem::EquationType::kDiffusion) {
    EXPECT_CALL(mock_builder, BuildDiffusionFormulation(Pointee(Ref(*finite_element_obs_ptr_)),
                                                        Pointee(Ref(*parameters.cross_sections_.value())),
                                                        _)).WillOnce(DoDefault());
    EXPECT_CALL(*diffusion_formulation_obs_ptr_, Precalculate(cells_.at(0)));
    using DiffusionFormulationPtr = std::unique_ptr<typename framework::builder::FrameworkBuilderI<dim>::DiffusionFormulation>;
    EXPECT_CALL(mock_builder, BuildUpdaterPointers(A<DiffusionFormulationPtr>(),
                                                   Pointee(Ref(*stamper_obs_ptr_)),
                                                   ContainerEq(reflective_boundaries))).WillOnce(DoDefault());
  }

  // End formulation specific calls, need_angular_storage should be set properly now

  EXPECT_CALL(mock_builder, BuildGroupSolution(n_angles)).WillOnce(DoDefault());

  EXPECT_CALL(mock_builder, BuildInitializer(Pointee(Ref(*updater_pointers_.fixed_updater_ptr)),
                                             parameters.neutron_energy_groups,
                                             n_angles)).WillOnce(DoDefault());
  EXPECT_CALL(*system_helper_mock_ptr_, SetUpMPIAngularSolution(Ref(*group_solution_obs_ptr_),
                                                                Ref(*domain_obs_ptr_),
                                                                1.0));
  EXPECT_CALL(mock_builder, BuildSingleGroupSolver(10000, 1e-10)).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildMomentConvergenceChecker(1e-6, 10000)).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildMomentMapConvergenceChecker(1e-6, 1000)).WillOnce(DoDefault());

  EXPECT_CALL(mock_builder, BuildGroupSolveIteration(Pointee(Ref(*single_group_solver_obs_ptr_)),
                                                     Pointee(Ref(*moment_convergence_checker_obs_ptr_)),
                                                     _,
                                                     Pointee(Ref(*group_solution_obs_ptr_)),
                                                     _,
                                                     Pointee(Ref(*moment_map_convergence_checker_obs_ptr_))))
      .WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildParameterConvergenceChecker(1e-6, 1000)).WillOnce(DoDefault());

  if (is_eigenvalue_solve) {
    if (parameters.k_effective_updater == eigenvalue::k_effective::K_EffectiveUpdaterName::kUpdaterViaRayleighQuotient) {
      EXPECT_CALL(mock_builder, BuildKEffectiveUpdater()).WillOnce(DoDefault());
    } else {
      EXPECT_CALL(mock_builder, BuildKEffectiveUpdater(Pointee(Ref(*finite_element_obs_ptr_)),
                                                       Pointee(Ref(*parameters.cross_sections_.value())),
                                                       Pointee(Ref(*domain_obs_ptr_)))).WillOnce(DoDefault());
    }
    EXPECT_CALL(mock_builder, BuildOuterIteration(Pointee(Ref(*group_solve_iteration_obs_ptr)),
                                                  Pointee(Ref(*parameter_convergence_checker_obs_ptr_)),
                                                  ::testing::NotNull(),
                                                  Pointee(Ref(*updater_pointers_.fission_source_updater_ptr)),
                                                  parameters.output_filename_base))
        .WillOnce(DoDefault());
  } else {
    EXPECT_CALL(mock_builder, BuildOuterIteration(Pointee(Ref(*group_solve_iteration_obs_ptr)),
                                                  Pointee(Ref(*parameter_convergence_checker_obs_ptr_)),
                                                  parameters.output_filename_base))
        .WillOnce(DoDefault());
  }

  std::set<FrameworkPart> needed_parts{ {FrameworkPart::ScatteringSourceUpdate} };
  if (is_eigenvalue_solve)
    needed_parts.insert(FrameworkPart::FissionSourceUpdate);
  if (need_angular_storage) {
    needed_parts.insert(FrameworkPart::AngularSolutionStorage);
    EXPECT_CALL(mock_validator_, AddPart(FrameworkPart::AngularSolutionStorage)).WillOnce(ReturnRef(mock_validator_));
  }
  EXPECT_CALL(mock_validator_, NeededParts()).WillOnce(Return(needed_parts));
  EXPECT_CALL(mock_validator_, ReportValidation());


  EXPECT_CALL(*group_solution_obs_ptr_, GetSolution(0)).WillOnce(DoDefault());
  EXPECT_CALL(mock_builder, BuildSystem(parameters.neutron_energy_groups,
                                        n_angles,
                                        Ref(*domain_obs_ptr_),
                                        _,
                                        is_eigenvalue_solve,
                                        need_angular_storage)).WillOnce(DoDefault());

  using ExpectedType = framework::Framework;
  auto framework_ptr = test_helper_ptr_->BuildFramework(this->mock_builder_, parameters);
  ASSERT_NE(framework_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<ExpectedType*>(framework_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  EXPECT_THAT(dynamic_ptr->system(), Pointee(Ref(*system_obs_ptr_)));
  if (is_eigenvalue_solve) {
    EXPECT_THAT(dynamic_ptr->outer_iterator_ptr(), Pointee(Ref(*outer_iteration_eigensolve_obs_ptr_)));
  } else {
    EXPECT_THAT(dynamic_ptr->outer_iterator_ptr(), Pointee(Ref(*outer_iteration_obs_ptr_)));
  }
  EXPECT_THAT(dynamic_ptr->initializer_ptr(), Pointee(Ref(*initializer_obs_ptr_)));
  EXPECT_THAT(dynamic_ptr->results_output_ptr(), NotNull());
}

TYPED_TEST_SUITE(FrameworkHelperBuildFrameworkIntegrationTests, bart::testing::AllDimensions);



// ===== BuildFramework ================================================================================================

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkDiffusion) {
  this->RunTest(this->default_parameters_);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkDiffusionEigensolve) {
  auto parameters{ this->default_parameters_ };
  parameters.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  this->RunTest(parameters);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAF) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
  this->RunTest(parameters);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAFEigensolve) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
  parameters.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  this->RunTest(parameters);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAFEigensolveRayleigh) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
  parameters.k_effective_updater = eigenvalue::k_effective::K_EffectiveUpdaterName::kUpdaterViaRayleighQuotient;
  parameters.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  this->RunTest(parameters);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAFReflectiveBCs) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
  parameters.reflective_boundaries = {problem::Boundary::kXMin};
  this->RunTest(parameters);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAFReflectiveBCsEigenSolve) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
  parameters.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  parameters.reflective_boundaries = {problem::Boundary::kXMin};
  this->RunTest(parameters);
}

} // namespace
