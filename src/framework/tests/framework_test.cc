#include "framework/framework.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class FrameworkTest : public ::testing::Test {
 public:
  using Framework = framework::Framework;


  std::unique_ptr<Framework> test_framework_;

  void SetUp() override;
};

void FrameworkTest::SetUp() {
  auto system_ptr = std::make_unique<system::System>();

  test_framework_ = std::make_unique<Framework>(
      std::move(system_ptr)
      );
}

TEST_F(FrameworkTest, Constructor) {
  EXPECT_NE(test_framework_->system(), nullptr);
}


} // namespace