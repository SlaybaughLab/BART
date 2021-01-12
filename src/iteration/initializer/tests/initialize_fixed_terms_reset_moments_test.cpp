#include "iteration/initializer/initialize_fixed_terms_reset_moments_test.hpp"

#include "formulation/updater/tests/fixed_updater_mock.h"
#include "system/moments/spherical_harmonic_types.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::ContainerEq, ::testing::DoDefault, ::testing::NiceMock, ::testing::Ref, ::testing::Return;

class InitializeFixedTermsResetMomentsTest : public ::testing::Test {
 public:
  using FixedUpdater = formulation::updater::FixedUpdaterMock;
  using TestInitializer = iteration::initializer::InitializeFixedTermsResetMoments;
  using SphericalHarmonics = NiceMock<system::moments::SphericalHarmonicMock>;

  // Test object
  std::unique_ptr<TestInitializer> test_initializer_ptr_{ nullptr };

  // Supporting objects
  system::System test_system_;

  // Supporting mocks
  std::shared_ptr<FixedUpdater> fixed_updater_ptr_{ nullptr };
  SphericalHarmonics* previous_moments_obs_ptr_, *current_moments_obs_ptr_;
  system::moments::MomentsMap previous_moments_map_, current_moments_map_;
  system::moments::MomentsMap reset_moments_map_;

  // Test parameters
  const int n_groups{ test_helpers::RandomInt(1, 5) };
  const int n_angles{ test_helpers::RandomInt(2, 3) };
  const int max_harmonic_l{ test_helpers::RandomInt(2, 3) };
  const int moment_vector_size{ test_helpers::RandomInt(10, 20) };

  // Test functions
  auto SetUp() -> void override;
};

auto InitializeFixedTermsResetMomentsTest::SetUp() -> void {
  for (int group = 0; group < n_groups; ++group) {
    for (int harmonic_l = 0; harmonic_l < max_harmonic_l; ++harmonic_l) {
      for (int harmonic_m = -harmonic_l; harmonic_m <= harmonic_l; ++harmonic_m) {
        dealii::Vector<double> previous_moment(moment_vector_size), current_moment(moment_vector_size);
        dealii::Vector<double> reset_moment(moment_vector_size);
        previous_moment = 2.0;
        current_moment = 2.0;
        reset_moment = 1.0;
        previous_moments_map_[{group, harmonic_l, harmonic_m}] = previous_moment;
        current_moments_map_[{group, harmonic_l, harmonic_m}] = current_moment;
        reset_moments_map_[{group, harmonic_l, harmonic_m}] = reset_moment;
      }
    }
  }
  test_system_.current_moments = std::make_unique<SphericalHarmonics>();
  test_system_.previous_moments = std::make_unique<SphericalHarmonics>();
  current_moments_obs_ptr_ = dynamic_cast<SphericalHarmonics*>(test_system_.current_moments.get());
  previous_moments_obs_ptr_ = dynamic_cast<SphericalHarmonics*>(test_system_.previous_moments.get());

  ON_CALL(*current_moments_obs_ptr_, begin()).WillByDefault(Return(current_moments_map_.begin()));
  ON_CALL(*current_moments_obs_ptr_, end()).WillByDefault(Return(current_moments_map_.end()));
  ON_CALL(*previous_moments_obs_ptr_, begin()).WillByDefault(Return(previous_moments_map_.begin()));
  ON_CALL(*previous_moments_obs_ptr_, end()).WillByDefault(Return(previous_moments_map_.end()));

  test_system_.k_effective = test_helpers::RandomDouble(0, 100);

  fixed_updater_ptr_ = std::make_shared<FixedUpdater>();
  test_initializer_ptr_ = std::make_unique<TestInitializer>(fixed_updater_ptr_, n_groups, n_angles);
}



TEST_F(InitializeFixedTermsResetMomentsTest, InitializeAndReset) {
  for (int group = 0; group < n_groups; ++group) {
    for (int angle = 0; angle < n_angles; ++angle) {
      system::EnergyGroup energy_group(group);
      quadrature::QuadraturePointIndex angle_index(angle);
      EXPECT_CALL(*fixed_updater_ptr_, UpdateFixedTerms(Ref(test_system_), energy_group, angle_index));
    }
  }
  EXPECT_CALL(*current_moments_obs_ptr_, begin()).WillOnce(DoDefault());
  EXPECT_CALL(*current_moments_obs_ptr_, end()).WillOnce(DoDefault());
  EXPECT_CALL(*previous_moments_obs_ptr_, begin()).WillOnce(DoDefault());
  EXPECT_CALL(*previous_moments_obs_ptr_, end()).WillOnce(DoDefault());

  test_initializer_ptr_->Initialize(test_system_);

  EXPECT_THAT(current_moments_map_, ContainerEq(reset_moments_map_));
  EXPECT_THAT(previous_moments_map_, ContainerEq(reset_moments_map_));
  EXPECT_EQ(test_system_.k_effective, std::nullopt);
}

} // namespace