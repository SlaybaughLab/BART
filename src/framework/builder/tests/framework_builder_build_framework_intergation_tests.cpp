#include "framework/builder/framework_builder.hpp"
#include "framework/builder/tests/framework_builder_mock.hpp"

#include "framework/framework_parameters.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class FrameworkBuilderBuildFrameworkIntegrationTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using FrameworkBuilder = framework::builder::FrameworkBuilder<dim>;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto FrameworkBuilderBuildFrameworkIntegrationTests<DimensionWrapper>::SetUp() -> void {}

TYPED_TEST_SUITE(FrameworkBuilderBuildFrameworkIntegrationTests, bart::testing::AllDimensions);

// ===== BuildFramework ================================================================================================

TYPED_TEST(FrameworkBuilderBuildFrameworkIntegrationTests, BuildFrameworkDefaultParameters) {
  EXPECT_TRUE(true);
}

} // namespace
