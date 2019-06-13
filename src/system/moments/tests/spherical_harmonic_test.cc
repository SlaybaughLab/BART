#include "system/moments/spherical_harmonic.h"

#include <array>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class SystemMomentsSphericalHarmonicTest : public ::testing::Test {
 protected:

  SystemMomentsSphericalHarmonicTest()
      : test_moments(total_groups, max_harmonic_l) {};

  // Test object
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
      EXPECT_ANY_THROW({
        system::moments::SphericalHarmonic bad_moment(group, max_harmonic_l);
      });
    }
  }
}


} // namespace