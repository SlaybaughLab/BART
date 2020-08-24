#include "utility/runtime/runtime_helper.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using ::testing::HasSubstr;

class UtilityRuntimeHelperTest : public ::testing::Test {
 public:
  const std::string version_{"7.6.9"};
  utility::runtime::RuntimeHelper test_helper{version_};
};

TEST_F(UtilityRuntimeHelperTest, Version) {
  EXPECT_EQ(test_helper.version(), version_);
}

TEST_F(UtilityRuntimeHelperTest, ProgramHeader) {
  auto header = test_helper.ProgramHeader();
  EXPECT_THAT(header, HasSubstr(version_));
}



} // namespace
