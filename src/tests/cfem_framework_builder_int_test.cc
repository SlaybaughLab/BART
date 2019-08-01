#include "framework/builder/cfem_framework_builder.h"
#include "problem/tests/parameters_mock.h"

#include <vector>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;

template <typename DimensionWrapper>
class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;

  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder<dim>;
  using ProblemParameters = problem::ParametersMock;

  FrameworkBuilder test_builder;
  ProblemParameters parameters;

  // Test Parameters
  const int polynomial_degree = 2;

  void SetUp() override;

};

template <typename DimensionWrapper>
void IntegrationTestCFEMFrameworkBuilder<DimensionWrapper>::SetUp() {
  std::vector<double> spatial_max;
  std::vector<int> n_cells;
  for (int i = 0; i < this->dim; ++i) {
    spatial_max.push_back(i + 10);
    n_cells.push_back(i + 10);
  }

  ON_CALL(parameters, NCells())
      .WillByDefault(Return(n_cells));
  ON_CALL(parameters, SpatialMax())
      .WillByDefault(Return(spatial_max));
  ON_CALL(parameters, FEPolynomialDegree())
      .WillByDefault(Return(polynomial_degree));

}

TYPED_TEST_CASE(IntegrationTestCFEMFrameworkBuilder,
                bart::testing::AllDimensions);

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildStamperTest) {

  EXPECT_CALL(this->parameters, NCells())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, SpatialMax())
      .WillOnce(DoDefault());
  EXPECT_CALL(this->parameters, FEPolynomialDegree())
      .WillOnce(DoDefault());

  auto stamper_ptr =
      this->test_builder.BuildStamper(&this->parameters, "1 2");
}

} // namespace