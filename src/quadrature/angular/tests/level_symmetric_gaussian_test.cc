#include "quadrature/angular/level_symmetric_gaussian.h"

#include <cmath>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class QuadratureAngularLevelSymmetricGaussianTest : public ::testing::Test {
 public:
  // Squared test function to verify integration properties, L = 1, m = 1
  double squared_integrated_function(quadrature::CartesianPosition<3> position) {
    const double x = position.get().at(0);
    const double y = position.get().at(1);
    const double radius = r(position);
    return 3.0/8.0/M_PI*(x*x + y*y)/(radius*radius);
  }
  // Radius calculation
  double r(quadrature::CartesianPosition<3> position) {
    auto &[x, y, z] = position.get();
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
  }
};

// Constructor should set order value properly.
TEST_F(QuadratureAngularLevelSymmetricGaussianTest, Constructor) {
  quadrature::angular::LevelSymmetricGaussian test_quadrature{quadrature::Order(4)};
  EXPECT_EQ(test_quadrature.order(), 4);
}

// Bad values of order (negative, odd, and too large) should throw an error.
TEST_F(QuadratureAngularLevelSymmetricGaussianTest, BadOrders) {
  std::array bad_orders{0, 1, 3, 5, 18, -1, -3, -18};

  for (auto bad_order : bad_orders) {
    EXPECT_ANY_THROW({
      quadrature::angular::LevelSymmetricGaussian test_quad{quadrature::Order(bad_order)};
    });
  }
}

// All orders should integrate a spherical harmonic exactly.
TEST_F(QuadratureAngularLevelSymmetricGaussianTest, Integration) {
  std::array orders{2, 4, 6, 8, 10, 12, 14, 16};
  for (const auto order : orders) {
    quadrature::angular::LevelSymmetricGaussian test_quad{quadrature::Order(order)};
    const auto quadrature_set = test_quad.GenerateSet();
    EXPECT_EQ(quadrature_set.size(), order*(order + 2)/8);

    double sum = 0;
    for (const auto& quadrature_point : quadrature_set) {
      auto& [position, weight] = quadrature_point;
      sum += weight.get()*squared_integrated_function(position);
    }
    EXPECT_NEAR(1.0, 8*sum, 1e-6);
  }
}

} // namespace