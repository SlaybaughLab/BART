#include "framework/builder/cfem_framework_builder.h"

#include <functional>
#include <numeric>
#include <vector>

#include <deal.II/fe/fe_q.h>
#include <solver/gmres.h>

#include "convergence/reporter/mpi_noisy.h"
#include "formulation/cfem_diffusion_stamper.h"
#include "data/cross_sections/material_cross_sections.hpp"
#include "material/tests/mock_material.h"
#include "problem/tests/parameters_mock.hpp"
#include "domain/domain.hpp"
#include "domain/finite_element/finite_element_gaussian.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "domain/tests/domain_mock.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/tests/angular_stamper_mock.h"
#include "formulation/cfem_saaf_stamper.h"
#include "iteration/initializer/set_fixed_terms_once.h"
#include "iteration/updater/source_updater_gauss_seidel.h"
#include "iteration/updater/angular_source_updater_gauss_seidel.h"
#include "iteration/updater/fixed_updater.h"
#include "iteration/updater/angular_fixed_updater.h"
#include "convergence/final_checker_or_n.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/moments/convergence_checker_l1_norm.hpp"
#include "solver/group/single_group_solver.h"
#include "quadrature/quadrature_set.hpp"
#include "quadrature/tests/quadrature_set_mock.hpp"
#include "test_helpers/dealii_test_domain.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;

template <typename DimensionWrapper>
class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder<dim>;
  using ProblemParameters = NiceMock<problem::ParametersMock>;
  using Material = NiceMock<btest::MockMaterial>;

  IntegrationTestCFEMFrameworkBuilder()
      : mock_material() {}


  FrameworkBuilder test_builder;
  ProblemParameters parameters;
  Material mock_material;

  // Test Parameters
  const int polynomial_degree = 2;
  std::vector<double> spatial_max;
  std::vector<int> n_cells;
  const int n_energy_groups = 3;
  std::array<int, 4> dofs_per_cell_by_dim_{1, 3, 9, 27};

  void SetUp() override;

};

template <typename DimensionWrapper>
void IntegrationTestCFEMFrameworkBuilder<DimensionWrapper>::SetUp() {
  for (int i = 0; i < this->dim; ++i) {
    spatial_max.push_back(10);
    n_cells.push_back(2);
  }

  std::map<problem::Boundary, bool> reflective_bcs{
      {problem::Boundary::kXMin, true},
      {problem::Boundary::kXMax, true},
      {problem::Boundary::kYMin, false},
      {problem::Boundary::kYMax, false},
      {problem::Boundary::kZMin, false},
      {problem::Boundary::kZMax, false},
  };

  ON_CALL(parameters, NEnergyGroups())
      .WillByDefault(Return(n_energy_groups));
  ON_CALL(parameters, NCells())
      .WillByDefault(Return(n_cells));
  ON_CALL(parameters, SpatialMax())
      .WillByDefault(Return(spatial_max));
  ON_CALL(parameters, FEPolynomialDegree())
      .WillByDefault(Return(polynomial_degree));
  ON_CALL(parameters, TransportModel())
      .WillByDefault(Return(problem::EquationType::kDiffusion));
  ON_CALL(parameters, ReflectiveBoundary())
      .WillByDefault(Return(reflective_bcs));
}

TYPED_TEST_CASE(IntegrationTestCFEMFrameworkBuilder,
                bart::testing::AllDimensions);

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildLSAngularQuadratureSet) {
  constexpr int dim = this->dim;
  const int order = 4;
  EXPECT_CALL(this->parameters, AngularQuad())
      .WillOnce(Return(problem::AngularQuadType::kLevelSymmetricGaussian));
  EXPECT_CALL(this->parameters, AngularQuadOrder())
      .WillOnce(Return(order));

  if (dim == 3) {

    auto quadrature_set = this->test_builder.BuildAngularQuadratureSet(
        &this->parameters);

    ASSERT_NE(nullptr, quadrature_set);
    using ExpectedType = quadrature::QuadratureSet<dim>;
    ASSERT_NE(nullptr,
              dynamic_cast<ExpectedType *>(quadrature_set.get()));
    EXPECT_EQ(quadrature_set->size(), order * (order + 2));
  } else {
    EXPECT_ANY_THROW({
         auto quadrature_set = this->test_builder.BuildAngularQuadratureSet(
             &this->parameters);
    });
  }
}



TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildConvergenceReporterTest) {
  using ExpectedType = convergence::reporter::MpiNoisy;
  
  auto convergence_reporter_ptr = this->test_builder.BuildConvergenceReporter();
  
  ASSERT_NE(nullptr, 
      dynamic_cast<ExpectedType*>(convergence_reporter_ptr.get()));
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildFiniteElementTest) {
  EXPECT_CALL(this->parameters, FEPolynomialDegree())
      .WillOnce(DoDefault());

  auto finite_element_ptr =
      this->test_builder.BuildFiniteElement(&this->parameters);

  // Test that correct type was returned
  auto gaussian_ptr = dynamic_cast<domain::finite_element::FiniteElementGaussian<this->dim>*>(
      finite_element_ptr.get());
  EXPECT_NE(gaussian_ptr, nullptr);

  // Test that correct values were used in instantiation
  EXPECT_EQ(finite_element_ptr->polynomial_degree(), this->polynomial_degree);
  auto dealii_finite_element_ptr = dynamic_cast<dealii::FE_Q<this->dim>*>(
      finite_element_ptr->finite_element());
  EXPECT_NE(dealii_finite_element_ptr, nullptr);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildDomainTest) {
  constexpr int dim = this->dim;
  auto finite_element_ptr =
      std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<dim>>>();

  EXPECT_CALL(this->parameters, NCells())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, SpatialMax())
      .WillOnce(DoDefault());

  auto test_domain_ptr =
      this->test_builder.BuildDomain(
          &this->parameters, finite_element_ptr, "1 1 2 2");

  EXPECT_NE(dynamic_cast<domain::Definition<this->dim>*>(test_domain_ptr.get()),
      nullptr);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildDiffusionStamper) {
  constexpr int dim = this->dim;

  auto finite_element_ptr = std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<dim>>>();
  auto cross_sections_ptr =
      std::make_shared<data::CrossSections>(this->mock_material);
  auto domain_ptr = std::make_shared<NiceMock<domain::DefinitionMock<dim>>>();

  EXPECT_CALL(this->parameters, TransportModel())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, ReflectiveBoundary())
      .WillOnce(DoDefault());

  auto test_stamper_ptr =
      this->test_builder.BuildStamper(
          &this->parameters,
          domain_ptr,
          finite_element_ptr,
          cross_sections_ptr);

  ASSERT_NE(test_stamper_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<formulation::CFEM_DiffusionStamper<this->dim>*>(
                test_stamper_ptr.get()));

  std::unordered_set<problem::Boundary> boundaries{problem::Boundary::kXMin,
                                                   problem::Boundary::kXMax};
  EXPECT_EQ(test_stamper_ptr->reflective_boundaries(), boundaries);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildAngularStamperBadEquation) {
  constexpr int dim = this->dim;

  auto finite_element_ptr = std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<dim>>>();
  auto cross_sections_ptr =
      std::make_shared<data::CrossSections>(this->mock_material);
  auto domain_ptr = std::make_shared<NiceMock<domain::DefinitionMock<dim>>>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  EXPECT_CALL(this->parameters, TransportModel())
      .WillOnce(Return(problem::EquationType::kDiffusion));

  EXPECT_ANY_THROW({
  this->test_builder.BuildAngularStamper(
      &this->parameters,
      domain_ptr,
      finite_element_ptr,
      cross_sections_ptr,
      quadrature_set_ptr);
  });
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder,
    BuildScalarFormulationBadEquationTypes) {
  using EquationType = problem::EquationType;
  std::vector<EquationType> bad_types{EquationType::kEvenParity,
                                      EquationType::kNone,
                                      EquationType::kSelfAdjointAngularFlux};

  constexpr int dim = this->dim;
  auto finite_element_ptr = std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<dim>>>();
  auto cross_sections_ptr =
      std::make_shared<data::CrossSections>(this->mock_material);
  auto domain_ptr = std::make_shared<NiceMock<domain::DefinitionMock<dim>>>();

  for (auto bad_type : bad_types) {
    EXPECT_CALL(this->parameters, TransportModel())
        .WillOnce(Return(bad_type));
    EXPECT_ANY_THROW({
      auto test_stamper_ptr =
          this->test_builder.BuildStamper(&this->parameters,
              domain_ptr, finite_element_ptr, cross_sections_ptr);});
  }
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildSingleGroupSolver) {
  using ExpectedType = solver::group::SingleGroupSolver;

  auto solver_ptr = this->test_builder.BuildSingleGroupSolver(100, 1e-12);

  ASSERT_NE(nullptr, solver_ptr);

  auto dynamic_ptr = dynamic_cast<ExpectedType*>(solver_ptr.get());
  ASSERT_NE(nullptr, dynamic_ptr);

  using ExpectedLinearSolverType = solver::GMRES;

  auto linear_solver_ptr = dynamic_cast<ExpectedLinearSolverType*>(
      dynamic_ptr->linear_solver_ptr());

  ASSERT_NE(nullptr, linear_solver_ptr);
  EXPECT_EQ(linear_solver_ptr->convergence_tolerance(), 1e-12);
  EXPECT_EQ(linear_solver_ptr->max_iterations(), 100);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildSourceUpdater) {

  auto stamper_ptr = std::make_shared<NiceMock<formulation::CFEM_StamperMock>>();

  auto test_source_updater =
      this->test_builder.BuildSourceUpdater(&this->parameters,
                                            stamper_ptr);
  using SourceUpdater =
      iteration::updater::SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;

  ASSERT_NE(test_source_updater, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<SourceUpdater*>(test_source_updater.get()));
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildAngularSourceUpdater) {
  constexpr int dim = this->dim;
  auto stamper_ptr = std::make_shared<NiceMock<formulation::AngularStamperMock<dim>>>();
  auto quadrature_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  auto test_source_updater =
      this->test_builder.BuildSourceUpdater(&this->parameters,
                                            stamper_ptr,
                                            quadrature_ptr);
  using SourceUpdater = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;

  ASSERT_NE(test_source_updater, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<SourceUpdater*>(test_source_updater.get()));
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildConvergenceChecker) {
  const double max_delta = 1e-4;
  const int max_iterations = 100;

  auto convergence_ptr =
      this->test_builder.BuildParameterConvergenceChecker(
          max_delta,
          max_iterations);

  using ParameterConvergenceChecker =
      convergence::FinalCheckerOrN<double, convergence::parameters::SingleParameterChecker>;


  ASSERT_NE(convergence_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<ParameterConvergenceChecker*>(convergence_ptr.get()));
  EXPECT_EQ(convergence_ptr->max_iterations(), max_iterations);

}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildMomentConvergenceChecker) {
  const double max_delta = 1e-4;
  const int max_iterations = 100;

  auto convergence_ptr =
      this->test_builder.BuildMomentConvergenceChecker(
          max_delta,
          max_iterations);

  using MomentConvergenceChecker =
  convergence::FinalCheckerOrN<system::moments::MomentVector,
                               convergence::moments::SingleMomentCheckerI>;


  ASSERT_NE(convergence_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<MomentConvergenceChecker*>(convergence_ptr.get()));
  EXPECT_EQ(convergence_ptr->max_iterations(), max_iterations);

}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildInitializer) {

  auto stamper_ptr = std::make_shared<formulation::CFEM_StamperMock>();

  EXPECT_CALL(this->parameters, NEnergyGroups())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, TransportModel())
      .WillOnce(DoDefault());
  auto iteration_ptr = this->test_builder.BuildInitializer(&this->parameters,
                                                           stamper_ptr);

  ASSERT_NE(nullptr, iteration_ptr);

  using ExpectedInitializerType = iteration::initializer::SetFixedTermsOnce;
  using ExpectedUpdaterType =
  iteration::updater::FixedUpdater<formulation::CFEMStamperI>;

  auto dynamic_iteration_ptr =
      dynamic_cast<ExpectedInitializerType*>(iteration_ptr.get());

  ASSERT_NE(nullptr, dynamic_iteration_ptr);
  ASSERT_NE(nullptr,
      dynamic_cast<ExpectedUpdaterType*>(dynamic_iteration_ptr->fixed_updater_ptr()));

  EXPECT_EQ(dynamic_iteration_ptr->total_groups(), this->n_energy_groups);
  EXPECT_EQ(dynamic_iteration_ptr->total_angles(), 1);
}

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildAngularInitializer) {
  constexpr int dim = this->dim;
  auto stamper_ptr = std::make_shared<formulation::AngularStamperMock<dim>>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();
  const int n_angles = 5;

  EXPECT_CALL(this->parameters, NEnergyGroups())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, TransportModel())
      .WillOnce(Return(problem::EquationType::kSelfAdjointAngularFlux));
  EXPECT_CALL(*quadrature_set_ptr, size())
      .WillOnce(Return(n_angles));

  auto iteration_ptr = this->test_builder.BuildInitializer(&this->parameters,
                                                           stamper_ptr,
                                                           quadrature_set_ptr);

  ASSERT_NE(nullptr, iteration_ptr);

  using ExpectedInitializerType = iteration::initializer::SetFixedTermsOnce;
  using ExpectedUpdaterType = iteration::updater::AngularFixedUpdater<formulation::AngularStamperI<dim>>;

  auto dynamic_iteration_ptr =
      dynamic_cast<ExpectedInitializerType*>(iteration_ptr.get());

  ASSERT_NE(nullptr, dynamic_iteration_ptr);
  ASSERT_NE(nullptr,
            dynamic_cast<ExpectedUpdaterType*>(dynamic_iteration_ptr->fixed_updater_ptr()));

  EXPECT_EQ(dynamic_iteration_ptr->total_groups(), this->n_energy_groups);
  EXPECT_EQ(dynamic_iteration_ptr->total_angles(), n_angles);
}

template <typename DimensionWrapper>
class IntegrationTestCFEMFrameworkBuilderDealiiDomain :
    public IntegrationTestCFEMFrameworkBuilder<DimensionWrapper>,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value>{
 public:
  void SetUp() override;
};

template <typename DimensionWrapper>
void IntegrationTestCFEMFrameworkBuilderDealiiDomain<DimensionWrapper>::SetUp() {
  IntegrationTestCFEMFrameworkBuilder<DimensionWrapper>::SetUp();
  this->SetUpDealii();
}

TYPED_TEST_SUITE(IntegrationTestCFEMFrameworkBuilderDealiiDomain,
                 bart::testing::AllDimensions);

TYPED_TEST(IntegrationTestCFEMFrameworkBuilderDealiiDomain, BuildSAAFStamper) {
  constexpr int dim = this->dim;

  auto finite_element_ptr = std::make_shared<NiceMock<domain::finite_element::FiniteElementMock<dim>>>();
  auto cross_sections_ptr =
      std::make_shared<data::CrossSections>(this->mock_material);
  auto domain_ptr = std::make_shared<NiceMock<domain::DefinitionMock<dim>>>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  EXPECT_CALL(this->parameters, TransportModel())
      .WillOnce(Return(problem::EquationType::kSelfAdjointAngularFlux));
  EXPECT_CALL(*domain_ptr, Cells())
      .WillOnce(Return(this->cells_));

  auto test_stamper_ptr =
      this->test_builder.BuildAngularStamper(
          &this->parameters,
          domain_ptr,
          finite_element_ptr,
          cross_sections_ptr,
          quadrature_set_ptr);

  ASSERT_NE(test_stamper_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<formulation::CFEM_SAAF_Stamper<this->dim>*>(
                test_stamper_ptr.get()));
}

} // namespace