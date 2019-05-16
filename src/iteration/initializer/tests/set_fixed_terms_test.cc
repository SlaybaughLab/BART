#include "iteration/initializer/set_fixed_terms.h"

#include "iteration/updater/tests/fixed_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::NiceMock;

class IterationInitializerSetFixedTermsTest : public ::testing::Test {
 protected:
  using InitializerType = iteration::initializer::SetFixedTerms;
  using MockFixedUpdaterType = iteration::updater::FixedUpdaterMock;

  std::unique_ptr<InitializerType> test_initializer_;
  MockFixedUpdaterType* updater_obs_ptr_ = nullptr;

  void SetUp() override;
};

void IterationInitializerSetFixedTermsTest::SetUp() {
  auto mock_fixed_updater_ptr =
      std::make_unique<NiceMock<MockFixedUpdaterType>>();

  test_initializer_ =
      std::make_unique<InitializerType>(std::move(mock_fixed_updater_ptr));

  updater_obs_ptr_ =
      dynamic_cast<MockFixedUpdaterType*>(test_initializer_->GetUpdater());
}

TEST_F(IterationInitializerSetFixedTermsTest, Constructor) {
  EXPECT_TRUE(updater_obs_ptr_ != nullptr);
}

} // namespace