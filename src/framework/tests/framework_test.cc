#include "framework/framework.h"

#include "iteration/initializer/tests/initializer_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class FrameworkTest : public ::testing::Test {
 public:
  using Framework = framework::Framework;
  using Initializer = iteration::initializer::InitializerMock;

  std::unique_ptr<Framework> test_framework_;

  void SetUp() override;
};

void FrameworkTest::SetUp() {
  auto system_ptr = std::make_unique<system::System>();
  auto initializer_ptr = std::make_unique<Initializer>();

  test_framework_ = std::make_unique<Framework>(
      std::move(system_ptr),
      std::move(initializer_ptr)
      );
}

TEST_F(FrameworkTest, Constructor) {
  EXPECT_NE(test_framework_->system(), nullptr);
  EXPECT_NE(test_framework_->initializer_ptr(), nullptr);
}


} // namespace