#include "quadrature/angular/level_symmetric_gaussian.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class QuadratureAngularLevelSymmetricGaussianTest : public ::testing::Test {};

TEST_F(QuadratureAngularLevelSymmetricGaussianTest, Constructor) {
  quadrature::angular::LevelSymmetricGaussian test_quadrature{quadrature::Order(4)};
  EXPECT_EQ(test_quadrature.order(), 4);
}

TEST_F(QuadratureAngularLevelSymmetricGaussianTest, BadOrders) {
  std::array bad_orders{0, 1, 3, 5, 18, -1, -3, -18};

  for (auto bad_order : bad_orders) {
    EXPECT_ANY_THROW({
      quadrature::angular::LevelSymmetricGaussian test_quad{quadrature::Order(bad_order)};
    });
  }
}

} // namespace