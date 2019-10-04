#include "quadrature/angular/level_symmetric_gaussian.h"

#include <cmath>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class QuadratureAngularLevelSymmetricGaussianTest : public ::testing::Test {
 public:

  double sh_y_0_1(quadrature::CartesianPosition<3> position) {
    double z = position.get().at(2);
    return 0.5*std::sqrt(3.0/M_PI)*z/r(position);
  }

  double r(quadrature::CartesianPosition<3> position) {
    auto &[x, y, z] = position.get();
    return std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
  }

};



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

TEST_F(QuadratureAngularLevelSymmetricGaussianTest, Integration) {
  std::array orders{2, 4, 6, 8, 10, 12, 14, 16};
  for (const auto order : orders) {
    quadrature::angular::LevelSymmetricGaussian test_quad{quadrature::Order(order)};

    const auto quadrature_set = test_quad.GenerateSet();
    EXPECT_EQ(quadrature_set.size(), order*(order + 2)/8);

    double sum = 0;
    for (const auto& quadrature_point : quadrature_set) {
      auto& [position, weight] = quadrature_point;
      sum += weight.get()*sh_y_0_1(position)*sh_y_0_1(position);
    }
    EXPECT_NEAR(1.0, 8*sum, 1e-6);

//    sum = 0;
//    for (const auto& quadrature_point : quadrature_set) {
//      auto& [position, weight] = quadrature_point;
//      sum += weight.get()*y_1_1_sq(position);
//    }
    //EXPECT_NEAR(1.0, 8*sum, 1e-6);
  }

}

} // namespace