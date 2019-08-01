#include "framework/builder/framework_builder.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IntegrationTestFrameworkBuilder : public ::testing::Test {
 protected:
  using FrameworkBuilder = framework::builder::FrameworkBuilder;

  FrameworkBuilder test_builder;
};

TEST_F(IntegrationTestFrameworkBuilder, Constructor) {
  EXPECT_TRUE(true);
}

} // namespace