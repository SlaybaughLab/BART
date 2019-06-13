#include "system/moments/spherical_harmonic.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class SystemMomentsSphericalHarmonicTest : public ::testing::Test {
 protected:

  // Test object
  system::moments::SphericalHarmonic test_moments;
};

TEST_F(SystemMomentsSphericalHarmonicTest, Dummy) {
  EXPECT_TRUE(true);
}


} // namespace