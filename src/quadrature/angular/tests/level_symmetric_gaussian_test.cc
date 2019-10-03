#include "quadrature/angular/level_symmetric_gaussian.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class QuadratureAngularLevelSymmetricGaussianTest : public ::testing::Test {};

TEST_F(QuadratureAngularLevelSymmetricGaussianTest, Dummy) {
  int order = 4;
  quadrature::angular::LevelSymmetricGaussian test_quadrature{quadrature::Order(order)};
  EXPECT_EQ(test_quadrature.order(), order);
}

} // namespace