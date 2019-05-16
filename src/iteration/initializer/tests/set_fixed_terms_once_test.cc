#include "iteration/initializer/set_fixed_terms_once.h"

#include <array>

#include "data/system.h"
#include "data/system/tests/bilinear_term_mock.h"
#include "iteration/updater/tests/fixed_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Ref, ::testing::Return;

class IterationInitializerSetFixedTermsOnceTest : public ::testing::Test {
 protected:
  using InitializerType = iteration::initializer::SetFixedTermsOnce;
  using MockFixedUpdaterType = iteration::updater::FixedUpdaterMock;
  using MockBilinearTermType = data::system::BilinearTermMock;

  IterationInitializerSetFixedTermsOnceTest()
      : total_groups_(btest::RandomDouble(1, 5)),
        total_angles_(btest::RandomDouble(1, 5)) {}

  // Initializer to be tested
  std::unique_ptr<InitializerType> test_initializer_;

  // Supporting objects
  data::System test_system_;

  // Pointers for observing mocks owned by other objects
  MockFixedUpdaterType* updater_obs_ptr_ = nullptr;

  // Random values for testing
  const int total_groups_;
  const int total_angles_;

  void SetUp() override;
};

void IterationInitializerSetFixedTermsOnceTest::SetUp() {
  // Set up testing object
  auto mock_fixed_updater_ptr =
      std::make_unique<MockFixedUpdaterType>();
  test_initializer_ =
      std::make_unique<InitializerType>(std::move(mock_fixed_updater_ptr),
                                        total_groups_,
                                        total_angles_);

  // Set up supporting objects
  auto mock_bilinear_term = std::make_unique<MockBilinearTermType>();
  test_system_.left_hand_side_ptr_ = std::move(mock_bilinear_term);

  // Set up observing pointers
  updater_obs_ptr_ =
      dynamic_cast<MockFixedUpdaterType*>(test_initializer_->fixed_updater_ptr());
}

TEST_F(IterationInitializerSetFixedTermsOnceTest, Constructor) {
  EXPECT_TRUE(updater_obs_ptr_ != nullptr);
  EXPECT_EQ(test_initializer_->total_angles(), total_angles_);
  EXPECT_EQ(test_initializer_->total_groups(), total_groups_);
}

TEST_F(IterationInitializerSetFixedTermsOnceTest, ConstructorThrows) {
  // Constructor should throw if fixed updater ptr is null
  std::unique_ptr<MockFixedUpdaterType> null_fixed_updater = nullptr;
  EXPECT_ANY_THROW(InitializerType initializer(std::move(null_fixed_updater),
                                               total_groups_, total_angles_););

  // Constructor should throw for bad groups and angle values
  std::array<int, 2> bad_values = {0, -1};

  for (int value : bad_values) {
    auto fixed_updater = std::make_unique<MockFixedUpdaterType>();
    EXPECT_ANY_THROW(InitializerType initializer(std::move(fixed_updater),
                                                 value, total_angles_););
    fixed_updater = std::make_unique<MockFixedUpdaterType>();
    EXPECT_ANY_THROW(InitializerType initializer(std::move(fixed_updater),
                                                 total_groups_, value););
  }
}

TEST_F(IterationInitializerSetFixedTermsOnceTest, Initialize) {
  // Initializer should access all left hand side terms (all groups/angles)
  for (int group = 0; group < total_groups_; ++group) {
    for (int angle = 0; angle < total_angles_; ++angle) {
      EXPECT_CALL(*updater_obs_ptr_, UpdateFixedTerms(Ref(test_system_),
                                                      group, angle));
    }
  }

  test_initializer_->Initialize(test_system_);
  test_initializer_->Initialize(test_system_);
}

} // namespace