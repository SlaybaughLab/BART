#include "iteration/initializer/initialize_fixed_terms.h"
#include "formulation/updater/tests/fixed_updater_mock.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class IterationInitializerInitializeFixedTermsTest : public ::testing::Test {
 public:
  using FixedUpdaterType = formulation::updater::FixedUpdaterMock;
  using InitializerType = iteration::initializer::InitializeFixedTerms;

  std::unique_ptr<InitializerType> test_initializer_ptr_;
  void SetUp() override;
};

void IterationInitializerInitializeFixedTermsTest::SetUp() {
  auto updater_ptr = std::make_unique<FixedUpdaterType>();
}

TEST_F(IterationInitializerInitializeFixedTermsTest, Constructor) {
  auto updater_ptr = std::make_unique<FixedUpdaterType>();
  EXPECT_NO_THROW({
    test_initializer_ptr_ = std::make_unique<InitializerType>(
        std::move(updater_ptr));
  });
  ASSERT_NE(test_initializer_ptr_->fixed_updater_ptr(), nullptr);
}

TEST_F(IterationInitializerInitializeFixedTermsTest, ConstructorBadDependency) {
  EXPECT_ANY_THROW({
  test_initializer_ptr_ = std::make_unique<InitializerType>(nullptr);
  });
}




} // namespace
