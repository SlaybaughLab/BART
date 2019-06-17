#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"

#include "system/solution/tests/mpi_angular_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class QuadCalcSphericalHarmonicMomentsOnlyScalar : public ::testing::Test {
 protected:
  // Tested object
  quadrature::calculators::SphericalHarmonicZerothMoment test_calculator;

  // Supporting objects
  system::solution::MPIAngularMock mock_solution;
};

TEST_F(QuadCalcSphericalHarmonicMomentsOnlyScalar, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace
