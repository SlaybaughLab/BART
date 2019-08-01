#include "framework/builder/cfem_framework_builder.h"
#include "problem/tests/parameters_mock.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;

  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder<dim>;
  using ProblemParameters = problem::ParametersMock;

  FrameworkBuilder test_builder;
  ProblemParameters parameters;

};

TYPED_TEST_CASE(IntegrationTestCFEMFrameworkBuilder,
                bart::testing::AllDimensions);

TYPED_TEST(IntegrationTestCFEMFrameworkBuilder, BuildStamperTest) {

  auto stamper_ptr =
      this->test_builder.BuildStamper(&this->parameters, "1 2");
}

} // namespace