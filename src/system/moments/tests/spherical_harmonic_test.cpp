#include "system/moments/spherical_harmonic.hpp"

#include <array>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Ref, ::testing::WhenDynamicCastTo, ::testing::NotNull;

class SystemMomentsSphericalHarmonicTest : public ::testing::Test {
 public:
  SystemMomentsSphericalHarmonicTest() : test_moments(total_groups, max_harmonic_l) {}

  // Test objects
  system::moments::SphericalHarmonic test_moments;

  // Test parameters
  static constexpr int max_harmonic_l = 2;
  static constexpr int total_groups = 2;
};

TEST_F(SystemMomentsSphericalHarmonicTest, Constructor) {
  EXPECT_EQ(test_moments.total_groups(), total_groups);
  EXPECT_EQ(test_moments.max_harmonic_l(), max_harmonic_l);
  int total_moments = total_groups * (max_harmonic_l + 1) * (max_harmonic_l + 1);
  EXPECT_EQ(test_moments.moments().size(), total_moments);
}

TEST_F(SystemMomentsSphericalHarmonicTest, BadGroupsAndHarmonics) {
  std::array<int, 2> bad_groups{0, -1};
  std::array<int, 1> bad_max_harmonic_l{-1};

  for (const auto group : bad_groups) {
    for (const auto max_harmonic_l : bad_max_harmonic_l) {
      EXPECT_ANY_THROW({ system::moments::SphericalHarmonic bad_moment(group, max_harmonic_l); });
    }
  }
}

TEST_F(SystemMomentsSphericalHarmonicTest, BracketOperator) {
  for (const auto& moment_pair : test_moments.moments()) {
    const auto& moment = test_moments[moment_pair.first];
    const auto& moment_from_get = test_moments.GetMoment(moment_pair.first);
    auto& non_const_moment_from_get = test_moments.GetMoment(moment_pair.first);
    EXPECT_THAT(moment, Ref(moment_pair.second));
    EXPECT_THAT(moment_from_get, Ref(moment_pair.second));
    EXPECT_THAT(non_const_moment_from_get, Ref(moment_pair.second));
  }

  const auto& const_test_moments = test_moments;

  for (const auto& moment_pair : const_test_moments.moments()) {
    const auto& moment = const_test_moments[moment_pair.first];
    const auto& const_moment = const_test_moments.GetMoment(moment_pair.first);
    EXPECT_THAT(moment, Ref(moment_pair.second));
    EXPECT_THAT(const_moment, Ref(moment_pair.second));
  }
}

TEST_F(SystemMomentsSphericalHarmonicTest, Assignment) {
  system::moments::MomentVector moment{10};
  system::moments::MomentIndex index{0,0,0};
  moment = 5;
  test_moments[index] = moment;
  EXPECT_EQ(test_moments[index], moment);
}

TEST_F(SystemMomentsSphericalHarmonicTest, BeginEndIterators) {
  EXPECT_EQ(test_moments.begin(), test_moments.moments().begin());
  EXPECT_EQ(test_moments.end(), test_moments.moments().end());
  EXPECT_EQ(test_moments.cbegin(), test_moments.moments().cbegin());
  EXPECT_EQ(test_moments.cend(), test_moments.moments().cend());
}

} // namespace