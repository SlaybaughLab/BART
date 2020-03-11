
#include <deal.II/fe/fe_q.h>

#include "framework/builder/framework_builder.h"

// Instantiated concerete classes
#include "convergence/reporter/mpi_noisy.h"
#include "domain/finite_element/finite_element_gaussian.h"
#include "quadrature/quadrature_set.h"

// Mock objects
#include "material/tests/mock_material.h"
#include "problem/tests/parameters_mock.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::Return, ::testing::NiceMock, ::testing::DoDefault;

template <typename DimensionWrapper>
class FrameworkBuilderIntegrationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FrameworkBuilder = framework::builder::FrameworkBuilder<dim>;
  using ProblemParameters = NiceMock<problem::ParametersMock>;
  using Material = NiceMock<btest::MockMaterial>;

  FrameworkBuilderIntegrationTest()
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
void FrameworkBuilderIntegrationTest<DimensionWrapper>::SetUp() {
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

TYPED_TEST_CASE(FrameworkBuilderIntegrationTest,
                bart::testing::AllDimensions);

TYPED_TEST(FrameworkBuilderIntegrationTest, BuildConvergenceReporterTest) {
  using ExpectedType = convergence::reporter::MpiNoisy;

  auto convergence_reporter_ptr = this->test_builder.BuildConvergenceReporter();

  ASSERT_NE(nullptr,
            dynamic_cast<ExpectedType*>(convergence_reporter_ptr.get()));
}

TYPED_TEST(FrameworkBuilderIntegrationTest, BuildFiniteElementTest) {
  constexpr int dim = this->dim;
  EXPECT_CALL(this->parameters, FEPolynomialDegree())
      .WillOnce(DoDefault());

  using ExpectedType = domain::finite_element::FiniteElementGaussian<dim>;

  auto finite_element_ptr = this->test_builder.BuildFiniteElement(this->parameters);
  auto gaussian_ptr = dynamic_cast<ExpectedType*>(finite_element_ptr.get());

  EXPECT_NE(gaussian_ptr, nullptr);
  EXPECT_EQ(finite_element_ptr->polynomial_degree(), this->polynomial_degree);
  auto dealii_finite_element_ptr = dynamic_cast<dealii::FE_Q<dim>*>(
      finite_element_ptr->finite_element());
  EXPECT_NE(dealii_finite_element_ptr, nullptr);
}

TYPED_TEST(FrameworkBuilderIntegrationTest, BuildLSAngularQuadratureSet) {
  constexpr int dim = this->dim;
  const int order = 4;
  EXPECT_CALL(this->parameters, AngularQuad())
      .WillOnce(Return(problem::AngularQuadType::kLevelSymmetricGaussian));
  EXPECT_CALL(this->parameters, AngularQuadOrder())
      .WillOnce(Return(order));

  if (dim == 3) {
    using ExpectedType = quadrature::QuadratureSet<dim>;
    auto quadrature_set = this->test_builder.BuildQuadratureSet(this->parameters);
    ASSERT_NE(nullptr, quadrature_set);
    ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(quadrature_set.get()));
    EXPECT_EQ(quadrature_set->size(), order * (order + 2));
  } else {
    EXPECT_ANY_THROW({
      auto quadrature_set = this->test_builder.BuildQuadratureSet(this->parameters);
    });
  }
}



} // namespace
