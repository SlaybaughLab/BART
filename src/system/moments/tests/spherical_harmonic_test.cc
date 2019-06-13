#include "system/moments/spherical_harmonic.h"

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
}


} // namespace