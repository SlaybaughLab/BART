#include "quadrature/angular/gauss_legendre.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

/* Constructor test for quadrature::angular::GaussLegendre class */
class QuadratureAngularGaussLegendreConstructorTest : public ::testing::Test {
 public:
  // Test parameters
  const int n_points = 4; // number of quadrature points
};

// Constructor given a valid number of points should not throw
TEST_F(QuadratureAngularGaussLegendreConstructorTest, Constructor) {
  EXPECT_NO_THROW(quadrature::angular::GaussLegendre test_quadrature(n_points));
}

// Constructor given an invalid number of points should throw
TEST_F(QuadratureAngularGaussLegendreConstructorTest, ConstructorBadNPoints) {
  auto bad_n_points = test_helpers::RandomVector(5, -5, 0);
  bad_n_points.push_back(0);
  for (const int n_points : bad_n_points) {
    EXPECT_ANY_THROW({
      quadrature::angular::GaussLegendre test_quadrature(n_points);
                     });
  }
}

/* Tests for the quadrature::angular::GaussLegendre class. */
class QuadratureAngularGaussLegendreTest
    : public QuadratureAngularGaussLegendreConstructorTest {
 public:
  QuadratureAngularGaussLegendreTest() : test_quadrature_(n_points) {}

  // Object to be tested
  quadrature::angular::GaussLegendre test_quadrature_;
};

// Order getter should return number of points
TEST_F(QuadratureAngularGaussLegendreTest, Order) {
  EXPECT_EQ(n_points, test_quadrature_.order());
}

/* Generate set should return a set of points and weights that exactly
 * integrate a spherical harmonic (Y_1_1) is used here */
TEST_F(QuadratureAngularGaussLegendreTest, GenerateSet) {
  auto set = test_quadrature_.GenerateSet();
  EXPECT_EQ(set.size(), n_points);

  // Generate negative x points to integrate over -1, 1
  for (auto point_pair : set) {
    point_pair.first.get().at(0) *= -1;
    set.push_back(point_pair);
  }

  double sum{0.0};

  // Y l = 1, m = 1
  auto legendre = [](const double point) {
    return 0.25 * 3 / M_PI * point * point;
  };

  for (const auto point_pair : set) {
    const auto point = point_pair.first.get().at(0);
    const auto weight = point_pair.second.get();
    sum += weight * legendre(point);
  }
  EXPECT_NEAR(1.0, sum, 1e-5);
}

} // namespace
