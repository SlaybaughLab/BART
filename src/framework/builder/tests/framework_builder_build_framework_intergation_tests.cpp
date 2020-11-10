#include "framework/builder/framework_builder_i.hpp"
#include "framework/builder/tests/framework_builder_mock.hpp"

#include "framework/framework_parameters.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class FrameworkBuilderBuildFrameworkIntegrationTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using FrameworkBuidler = framework::builder::FrameworkBuilderMock<dim>;
  using FrameworkParameters = framework::FrameworkParameters;
  FrameworkBuidler mock_builder_;
  FrameworkParameters parameters_;
  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto FrameworkBuilderBuildFrameworkIntegrationTests<DimensionWrapper>::SetUp() -> void {}

TYPED_TEST_SUITE(FrameworkBuilderBuildFrameworkIntegrationTests, bart::testing::AllDimensions);

// ===== BuildFramework ================================================================================================

TYPED_TEST(FrameworkBuilderBuildFrameworkIntegrationTests, BuildFrameworkDefaultParameters) {
  auto framework_ptr = framework::builder::BuildFramework(this->mock_builder_, this->parameters_);
  ASSERT_NE(framework_ptr, nullptr);
}

} // namespace
