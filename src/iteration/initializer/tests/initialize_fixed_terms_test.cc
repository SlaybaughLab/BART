#include "iteration/initializer/initialize_fixed_terms.h"
#include "formulation/updater/tests/fixed_updater_mock.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

using namespace bart;

class IterationInitializerInitializeFixedTermsTest : public ::testing::Test {
 public:
  using FixedUpdaterType = formulation::updater::FixedUpdaterMock;
  using InitializerType = iteration::initializer::InitializeFixedTerms;

  std::unique_ptr<InitializerType> test_initializer_ptr_;

  // Test parameters
  const int total_groups_ = test_helpers::RandomDouble(1, 5);
  const int total_angles_ = test_helpers::RandomDouble(1, 5);

  void SetUp() override;
};

void IterationInitializerInitializeFixedTermsTest::SetUp() {
  auto updater_ptr = std::make_unique<FixedUpdaterType>();
}

TEST_F(IterationInitializerInitializeFixedTermsTest, Constructor) {
  auto updater_ptr = std::make_unique<FixedUpdaterType>();
  EXPECT_NO_THROW({
    test_initializer_ptr_ = std::make_unique<InitializerType>(
        std::move(updater_ptr), total_groups_, total_angles_);
  });
  ASSERT_NE(test_initializer_ptr_->fixed_updater_ptr(), nullptr);
  EXPECT_EQ(total_groups_, test_initializer_ptr_->total_groups());
  EXPECT_EQ(total_angles_, test_initializer_ptr_->total_angles());
}

TEST_F(IterationInitializerInitializeFixedTermsTest, ConstructorBadDependency) {
  EXPECT_ANY_THROW({
    test_initializer_ptr_ = std::make_unique<InitializerType>(nullptr,
                                                              total_groups_,
                                                              total_angles_);
  });
}

TEST_F(IterationInitializerInitializeFixedTermsTest, ConstructorBadGroupAngle) {
  std::array<int, 2> bad_values{0, -1};
  for (const int bad_value : bad_values) {
    EXPECT_ANY_THROW({
      auto updater_ptr = std::make_unique<FixedUpdaterType>();
      test_initializer_ptr_ = std::make_unique<InitializerType>(
          std::move(updater_ptr), total_groups_, bad_value);
    });
    EXPECT_ANY_THROW({
      auto updater_ptr = std::make_unique<FixedUpdaterType>();
      test_initializer_ptr_ = std::make_unique<InitializerType>(
          std::move(updater_ptr), bad_value, total_angles_);
    });
    EXPECT_ANY_THROW({
      auto updater_ptr = std::make_unique<FixedUpdaterType>();
      test_initializer_ptr_ = std::make_unique<InitializerType>(
          std::move(updater_ptr), bad_value, bad_value);
    });
  }
}




} // namespace
