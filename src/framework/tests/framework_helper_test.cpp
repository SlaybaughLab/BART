#include "framework/framework_helper.hpp"

#include "formulation/tests/stamper_mock.h"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/updater/tests/boundary_conditions_updater_mock.h"
#include "formulation/updater/tests/fission_source_updater_mock.h"
#include "formulation/updater/tests/fixed_updater_mock.h"
#include "formulation/updater/tests/scattering_source_updater_mock.h"
#include "formulation/scalar/tests/diffusion_mock.h"
#include "framework/builder/framework_builder_i.hpp"
#include "framework/builder/tests/framework_builder_mock.hpp"
#include "framework/framework_parameters.hpp"
#include "iteration/initializer/tests/initializer_mock.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "domain/tests/definition_mock.h"
#include "material/tests/mock_material.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "quadrature/calculators/tests/spherical_harmonic_moments_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "system/tests/system_helper_mock.hpp"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"

namespace  {

using namespace bart;
using ::testing::Return, ::testing::ByMove, ::testing::DoDefault, ::testing::_, ::testing::NiceMock;
using ::testing::Ref, ::testing::Pointee, ::testing::ReturnRef, ::testing::ContainerEq, ::testing::SizeIs;
using ::testing::A, ::testing::AllOf;
using ::testing::NotNull, ::testing::WhenDynamicCastTo;

template <typename DimensionWrapper>
class FrameworkHelperConstructorTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(FrameworkHelperConstructorTests, bart::testing::AllDimensions);

TYPED_TEST(FrameworkHelperConstructorTests, Constructor) {
  constexpr int dim = this->dim;
  using SystemHelper = const system::SystemHelperMock<dim>;
  auto system_helper_ptr = std::make_shared<SystemHelper>();
  framework::FrameworkHelper<dim> helper(system_helper_ptr);
  EXPECT_THAT(helper.system_helper_ptr(), WhenDynamicCastTo<SystemHelper*>(NotNull()));
}

TYPED_TEST(FrameworkHelperConstructorTests, ConstructorNullDependency) {
  constexpr int dim = this->dim;
  EXPECT_ANY_THROW({
    framework::FrameworkHelper<dim> helper(nullptr);
  });

}

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
  using InitializerMock = iteration::initializer::InitializerMock;
  using MomentCalculatorMock = quadrature::calculators::SphericalHarmonicMomentsMock;
  using QuadratureSetMock = typename quadrature::QuadratureSetMock<dim>;
  using StamperMock = typename formulation::StamperMock<dim>;
  using SAAFFormulationMock = typename formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using SystemHelper = typename system::SystemHelperMock<dim>;

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
  InitializerMock* initializer_obs_ptr_{ nullptr };
  MomentCalculatorMock* moment_calculator_obs_ptr_{ nullptr };
  std::shared_ptr<QuadratureSetMock> quadrature_set_mock_ptr_{ nullptr };
  SAAFFormulationMock* saaf_formulation_obs_ptr_{ nullptr };
  StamperMock* stamper_obs_ptr_{ nullptr };
  std::shared_ptr<SystemHelper> system_helper_mock_ptr_{ nullptr };

  UpdaterPointers updater_pointers_;

  // Test parameters and supporting objects
  FrameworkBuidler mock_builder_;
  FrameworkParameters default_parameters_;
  const int total_quadrature_angles{ test_helpers::RandomInt(10, 20) };
  std::vector<domain::CellPtr<dim>> cells_;

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
  auto initializer_ptr = std::make_unique<NiceMock<InitializerMock>>();
  initializer_obs_ptr_ = initializer_ptr.get();
  auto moment_calculator_ptr = std::make_unique<NiceMock<MomentCalculatorMock>>();
  moment_calculator_obs_ptr_ = moment_calculator_ptr.get();
  quadrature_set_mock_ptr_ = std::make_shared<QuadratureSetMock>();
  auto saaf_ptr = std::make_unique<NiceMock<SAAFFormulationMock>>();
  saaf_formulation_obs_ptr_ = saaf_ptr.get();
  auto stamper_ptr = std::make_unique<NiceMock<StamperMock>>();
  stamper_obs_ptr_ = stamper_ptr.get();
  system_helper_mock_ptr_ = std::make_shared<NiceMock<SystemHelper>>();

  updater_pointers_.boundary_conditions_updater_ptr = std::make_shared<NiceMock<BoundaryConditionsUpdaterMock>>();
  updater_pointers_.fission_source_updater_ptr = std::make_shared<NiceMock<FissionSourceUpdaterMock>>();
  updater_pointers_.fixed_updater_ptr = std::make_shared<NiceMock<FixedTermUpdaterMock>>();
  updater_pointers_.scattering_source_updater_ptr = std::make_shared<NiceMock<ScatteringSourceUpdaterMock>>();

  using DiffusionFormulationPtr = std::unique_ptr<typename FrameworkBuidler::DiffusionFormulation>;
  using SAAFFormulationPtr = std::unique_ptr<typename FrameworkBuidler::SAAFFormulation>;

  ON_CALL(mock_builder_, BuildDiffusionFormulation(_,_,_)).WillByDefault(ReturnByMove(diffusion_formulation_ptr));
  ON_CALL(mock_builder_, BuildDomain(_, _, _, _)).WillByDefault(ReturnByMove(domain_ptr));
  ON_CALL(mock_builder_, BuildFiniteElement(_,_,_)).WillByDefault(ReturnByMove(finite_element_ptr));
  ON_CALL(mock_builder_, BuildGroupSolution(_)).WillByDefault(ReturnByMove(group_solution_ptr));
  ON_CALL(mock_builder_, BuildInitializer(_,_,_)).WillByDefault(ReturnByMove(initializer_ptr));
  ON_CALL(mock_builder_, BuildQuadratureSet(_,_)).WillByDefault(Return(quadrature_set_mock_ptr_));
  ON_CALL(mock_builder_, BuildSAAFFormulation(_,_,_,_)).WillByDefault(ReturnByMove(saaf_ptr));
  ON_CALL(mock_builder_, BuildStamper(_)).WillByDefault(ReturnByMove(stamper_ptr));
  ON_CALL(mock_builder_, BuildUpdaterPointers(A<SAAFFormulationPtr>(),_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildUpdaterPointers(A<DiffusionFormulationPtr>(),_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildUpdaterPointers(_,_,_,_,_)).WillByDefault(Return(updater_pointers_));
  ON_CALL(mock_builder_, BuildMomentCalculator(_)).WillByDefault(ReturnByMove(moment_calculator_ptr));
  ON_CALL(mock_builder_, BuildMomentCalculator(_,_)).WillByDefault(ReturnByMove(moment_calculator_ptr));

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

  // Mock Builder calls
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
      EXPECT_CALL(*system_helper_mock_ptr_, SetUpEnergyGroupToAngularSolutionPtrMap(
          _, parameters.neutron_energy_groups, total_quadrature_angles));
      EXPECT_CALL(mock_builder, BuildUpdaterPointers(
          Pointee(Ref(*saaf_formulation_obs_ptr_)),
          Pointee(Ref(*stamper_obs_ptr_)),
          Pointee(Ref(*quadrature_set_mock_ptr_)),
          ContainerEq(reflective_boundaries),
          _))
          .WillOnce(DoDefault());
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

  EXPECT_CALL(mock_builder, BuildGroupSolution(n_angles)).WillOnce(DoDefault());

  EXPECT_CALL(mock_builder, BuildInitializer(Pointee(Ref(*updater_pointers_.fixed_updater_ptr)),
                                             parameters.neutron_energy_groups,
                                             n_angles)).WillOnce(DoDefault());


  auto framework_ptr = test_helper_ptr_->BuildFramework(this->mock_builder_, parameters);
  ASSERT_NE(framework_ptr, nullptr);
}

TYPED_TEST_SUITE(FrameworkHelperBuildFrameworkIntegrationTests, bart::testing::AllDimensions);



// ===== BuildFramework ================================================================================================

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkDefaultParameters) {
  this->RunTest(this->default_parameters_);
}

TYPED_TEST(FrameworkHelperBuildFrameworkIntegrationTests, BuildFrameworkSAAF) {
  auto parameters{ this-> default_parameters_ };
  using Order = framework::FrameworkParameters::AngularQuadratureOrder;
  parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  parameters.angular_quadrature_type = problem::AngularQuadType::kLevelSymmetricGaussian;
  parameters.angular_quadrature_order = Order(test_helpers::RandomInt(5, 10));
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

} // namespace
