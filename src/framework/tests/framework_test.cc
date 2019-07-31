#include "framework/framework.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class FrameworkTest : public ::testing::Test {
 public:

  void SetUp() override;
};

void FrameworkTest::SetUp() {

}

TEST_F(FrameworkTest, Constructor) {
  EXPECT_NO_THROW({
    framework::Framework new_framework;
  });
}


} // namespace