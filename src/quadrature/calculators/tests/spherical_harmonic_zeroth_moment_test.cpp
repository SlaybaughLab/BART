#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class QuadCalcSphericalHarmonicMomentsOnlyScalar : public ::testing::Test {
 protected:
  quadrature::calculators::SphericalHarmonicZerothMoment test_calculator;
};

TEST_F(QuadCalcSphericalHarmonicMomentsOnlyScalar, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace
