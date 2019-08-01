#include "framework/builder/cfem_framework_builder.h"
#include "problem/tests/parameters_mock.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IntegrationTestCFEMFrameworkBuilder : public ::testing::Test {
 protected:
  using FrameworkBuilder = framework::builder::CFEM_FrameworkBuilder;
  using ProblemParameters = problem::ParametersMock;

  FrameworkBuilder test_builder;
  ProblemParameters parameters;

};

TEST_F(IntegrationTestCFEMFrameworkBuilder, BuildStamperTest) {
  auto stamper_ptr = test_builder.BuildStamper(&parameters);
}

} // namespace