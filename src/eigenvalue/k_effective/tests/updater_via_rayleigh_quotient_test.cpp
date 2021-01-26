#include "eigenvalue/k_effective/updater_via_rayleigh_quotient.hpp"

#include "system/system.h"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::ReturnRef;

class K_EffectiveUpdaterViaRayleighQuotientTest : public ::testing::Test {
 public:
  using SphericalHarmonicMock = NiceMock<system::moments::SphericalHarmonicMock>;
  using Vector = dealii::Vector<double>;

  static constexpr int total_groups_{ 2 };
  static constexpr int vector_size_{ 3 };
  const double initial_k_effective_{ 1.0235 };
  const double expected_k_effective_{ 0.850701299 };
  std::array<Vector, total_groups_> current_moments_, previous_moments_;

  system::System test_system_ { .k_effective = initial_k_effective_, .total_groups = total_groups_ };

  SphericalHarmonicMock* current_moments_obs_ptr_;
  SphericalHarmonicMock* previous_moments_obs_ptr_;

  auto SetUp() -> void override;
};

auto K_EffectiveUpdaterViaRayleighQuotientTest::SetUp() -> void {
  for (int group = 0; group < total_groups_; ++group) {
    Vector current_group_moment(vector_size_), previous_group_moment(vector_size_);
    for (int i = 0; i < vector_size_; ++i) {
      if (group == 0) {
        current_group_moment[i] = (i + 1);
        previous_group_moment[i] = (i + 1 + vector_size_);
      } else {
        current_group_moment[i] = vector_size_ - i;
        previous_group_moment[i] = 2 * vector_size_ - i;
      }
    }
    current_moments_.at(group) = current_group_moment;
    previous_moments_.at(group) = previous_group_moment;
  }

  test_system_.current_moments = std::make_unique<SphericalHarmonicMock>();
  test_system_.previous_moments = std::make_unique<SphericalHarmonicMock>();
  current_moments_obs_ptr_ = dynamic_cast<SphericalHarmonicMock*>(test_system_.current_moments.get());
  previous_moments_obs_ptr_ = dynamic_cast<SphericalHarmonicMock*>(test_system_.previous_moments.get());

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    ON_CALL(*current_moments_obs_ptr_, GetMoment(index)).WillByDefault(ReturnRef(current_moments_.at(group)));
    ON_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillByDefault(ReturnRef(previous_moments_.at(group)));
  }
}

TEST_F(K_EffectiveUpdaterViaRayleighQuotientTest, Calculate) {
  eigenvalue::k_effective::UpdaterViaRayleighQuotient test_updater;
  EXPECT_FALSE(test_updater.k_effective().has_value());

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    EXPECT_CALL(*current_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
    EXPECT_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
  }

  auto calculated_k_effective = test_updater.CalculateK_Effective(test_system_);
  EXPECT_NEAR(calculated_k_effective, expected_k_effective_, 1e-8);
  ASSERT_TRUE(test_updater.k_effective().has_value());
  EXPECT_EQ(calculated_k_effective, test_updater.k_effective().value());
}

TEST_F(K_EffectiveUpdaterViaRayleighQuotientTest, ZeroCurrentFlux) {
  eigenvalue::k_effective::UpdaterViaRayleighQuotient test_updater;
  EXPECT_FALSE(test_updater.k_effective().has_value());

  Vector zero_flux(vector_size_);

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    EXPECT_CALL(*current_moments_obs_ptr_, GetMoment(index)).WillOnce(ReturnRef(zero_flux));
    EXPECT_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
  }

  auto calculated_k_effective = test_updater.CalculateK_Effective(test_system_);
  EXPECT_NEAR(calculated_k_effective, 0, 1e-8);
  ASSERT_TRUE(test_updater.k_effective().has_value());
  EXPECT_EQ(calculated_k_effective, test_updater.k_effective().value());
}

TEST_F(K_EffectiveUpdaterViaRayleighQuotientTest, ZeroPreviousFlux) {
  eigenvalue::k_effective::UpdaterViaRayleighQuotient test_updater;
  EXPECT_FALSE(test_updater.k_effective().has_value());

  Vector zero_flux(vector_size_);

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    EXPECT_CALL(*current_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
    EXPECT_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillOnce(ReturnRef(zero_flux));
  }

  auto calculated_k_effective = test_updater.CalculateK_Effective(test_system_);
  EXPECT_NEAR(calculated_k_effective, initial_k_effective_, 1e-8);
  ASSERT_TRUE(test_updater.k_effective().has_value());
  EXPECT_EQ(calculated_k_effective, test_updater.k_effective().value());
}

TEST_F(K_EffectiveUpdaterViaRayleighQuotientTest, ZeroPreviousFluxOneGroup) {
  eigenvalue::k_effective::UpdaterViaRayleighQuotient test_updater;
  EXPECT_FALSE(test_updater.k_effective().has_value());

  Vector zero_flux(vector_size_);

  for (int group = 0; group < total_groups_; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    EXPECT_CALL(*current_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
    if (group == 0) {
      EXPECT_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillOnce(ReturnRef(zero_flux));
    } else {
      EXPECT_CALL(*previous_moments_obs_ptr_, GetMoment(index)).WillOnce(DoDefault());
    }
  }

  auto calculated_k_effective = test_updater.CalculateK_Effective(test_system_);
  EXPECT_NEAR(calculated_k_effective, 0.5*expected_k_effective_, 1e-8);
  ASSERT_TRUE(test_updater.k_effective().has_value());
  EXPECT_EQ(calculated_k_effective, test_updater.k_effective().value());
}

} // namespace
