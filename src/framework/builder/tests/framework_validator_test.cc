#include "test_helpers/gmock_wrapper.h"

#include "framework/builder/framework_validator.h"
#include "problem/tests/parameters_mock.h"
#include "utility/reporter/tests/basic_reporter_mock.h"

namespace  {

using namespace bart;

class FrameworkBuilderFrameworkValidatorTest : public ::testing::Test {
 public:
  framework::builder::FrameworkValidator test_validator;
};

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTest) {

}

} // namespace
