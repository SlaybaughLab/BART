#include "framework/builder/cfem_framework_builder.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 protected:
  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder;

  FrameworkBuilder test_builder;
};

TEST_F(IntegrationTestCFEMFrameworkBuilder, Constructor) {
  EXPECT_TRUE(true);
}

} // namespace