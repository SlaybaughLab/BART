#include "iteration/initializer/initialize_fixed_terms_once.h"

#include <array>

#include "system/system.hpp"
#include "formulation/updater/tests/fixed_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::NiceMock, ::testing::Ref, ::testing::Return;

class IterationInitializerInitializeFixedTermsOnceTest : public ::testing::Test {
 protected:
  using InitializerType = iteration::initializer::InitializeFixedTermsOnce;
  using MockFixedUpdaterType = formulation::updater::FixedUpdaterMock;

  IterationInitializerInitializeFixedTermsOnceTest()
      : total_groups_(test_helpers::RandomDouble(1, 5)),
        total_angles_(test_helpers::RandomDouble(1, 5)) {}

  // Initializer to be tested
  std::unique_ptr<InitializerType> test_initializer_;

  // Supporting objects
  system::System test_system_;

  // Pointers for observing mocks owned by other objects
  MockFixedUpdaterType* updater_obs_ptr_ = nullptr;

  // Random values for testing
  const int total_groups_;
  const int total_angles_;

  void SetUp() override;
};

void IterationInitializerInitializeFixedTermsOnceTest::SetUp() {
  // Set up testing object
  auto mock_fixed_updater_ptr =
      std::make_shared<MockFixedUpdaterType>();
  updater_obs_ptr_ = mock_fixed_updater_ptr.get();

  test_initializer_ =
      std::make_unique<InitializerType>(mock_fixed_updater_ptr,
                                        total_groups_,
                                        total_angles_);
}

TEST_F(IterationInitializerInitializeFixedTermsOnceTest, Constructor) {
  EXPECT_TRUE(updater_obs_ptr_ != nullptr);
  EXPECT_EQ(test_initializer_->total_angles(), total_angles_);
  EXPECT_EQ(test_initializer_->total_groups(), total_groups_);
}

TEST_F(IterationInitializerInitializeFixedTermsOnceTest, ConstructorThrows) {
  // Constructor should throw if fixed updater ptr is null
  std::shared_ptr<MockFixedUpdaterType> null_fixed_updater = nullptr;
  EXPECT_ANY_THROW(InitializerType initializer(null_fixed_updater,
                                               total_groups_, total_angles_););

  // Constructor should throw for bad groups and angle values
  std::array<int, 2> bad_values = {0, -1};

  for (int value : bad_values) {
    auto fixed_updater = std::make_shared<MockFixedUpdaterType>();
    EXPECT_ANY_THROW(InitializerType initializer(fixed_updater,
                                                 value, total_angles_););
    EXPECT_ANY_THROW(InitializerType initializer(fixed_updater,
                                                 total_groups_, value););
  }
}

TEST_F(IterationInitializerInitializeFixedTermsOnceTest, Initialize) {
  // Initializer should access all left hand side terms (all groups/angles)
  // This will run only twice despite us calling it three times.
  for (int group = 0; group < total_groups_; ++group) {
    for (int angle = 0; angle < total_angles_; ++angle) {
      system::EnergyGroup energy_group(group);
      quadrature::QuadraturePointIndex angle_index(angle);
      EXPECT_CALL(*updater_obs_ptr_, UpdateFixedTerms(Ref(test_system_),
                                                      energy_group,
                                                      angle_index))
          .Times(2);
    }
  }
  EXPECT_FALSE(test_initializer_->initialize_was_called());
  test_initializer_->Initialize(test_system_);

  // Initialize shouldn't do anything because initialize has been called
  EXPECT_TRUE(test_initializer_->initialize_was_called());
  test_initializer_->Initialize(test_system_);

  // Reset status of initialize called
  test_initializer_->set_initialize_was_called(false);
  test_initializer_->Initialize(test_system_);
  EXPECT_TRUE(test_initializer_->initialize_was_called());
}

} // namespace