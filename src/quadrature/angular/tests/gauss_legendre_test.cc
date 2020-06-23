#include "quadrature/angular/gauss_legendre.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class QuadratureAngularGaussLegendreTest : public ::testing::Test {
 public:
  QuadratureAngularGaussLegendreTest() : test_quadrature_(n_points) {}
  quadrature::angular::GaussLegendre test_quadrature_;
  const int n_points = 4;

};

TEST_F(QuadratureAngularGaussLegendreTest, Constructor) {
  EXPECT_NO_THROW(quadrature::angular::GaussLegendre test_quadrature(n_points));
}

TEST_F(QuadratureAngularGaussLegendreTest, ConstructorBadNPoints) {
  auto bad_n_points = test_helpers::RandomVector(5, -5, 0);
  bad_n_points.push_back(0);
  for (const int n_points : bad_n_points) {
    EXPECT_ANY_THROW({
      quadrature::angular::GaussLegendre test_quadrature(n_points);
                     });
  }
}

TEST_F(QuadratureAngularGaussLegendreTest, Order) {
  EXPECT_EQ(n_points, test_quadrature_.order());
}

TEST_F(QuadratureAngularGaussLegendreTest, GenerateSet) {
  auto set = test_quadrature_.GenerateSet();
  EXPECT_EQ(set.size(), n_points);

  double sum{0.0};

  // P_1 Legendre polynomial should integrate to 0.4 over [-1, 1]
  auto legendre = [](const double point) {
    return 0.5 * (3*point*point - 1);
  };

  for (const auto point_pair : set) {
    const auto point = point_pair.first.get().at(0);
    const auto weight = point_pair.second.get();
    sum += weight * legendre(point) * legendre(point);
  }
  // Multiply by 2 because points are generated from [0, 1]
  EXPECT_NEAR(0.4, 2*sum, 1e-5);
}

} // namespace
