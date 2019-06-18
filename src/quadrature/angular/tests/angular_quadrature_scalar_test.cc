#include "quadrature/angular/angular_quadrature_scalar.h"

#include <array>
#include <map>
#include <vector>

#include "system/system_types.h"
#include "quadrature/angular/angular_quadrature_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::ContainerEq;

template <typename DimensionWrapper>
class AngularQuadratureScalarTests : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  // Tested object
  quadrature::angular::AngularQuadratureScalar<dim> test_scalar_quad_;
};

TYPED_TEST_CASE(AngularQuadratureScalarTests, bart::testing::AllDimensions);

TYPED_TEST(AngularQuadratureScalarTests, Constructor) {
  auto& test_scalar_quad_ = this->test_scalar_quad_;
  constexpr int dim = this->dim;
  using QuadraturePoint = quadrature::angular::QuadraturePoint<dim>;
  using Ordinate = quadrature::angular::Ordinate<dim>;
  using Weight = quadrature::angular::Weight;
  using system::AngleIndex;

  EXPECT_EQ(test_scalar_quad_.total_quadrature_points(), 1);

  Ordinate zero_ordinate;

  for (auto& angle_coordinate : zero_ordinate)
    angle_coordinate = 0;

  std::map<AngleIndex, QuadraturePoint> expected_quad_pt_map{{0, {1.0, zero_ordinate}}};
  std::vector<QuadraturePoint> expected_quadrature_points{{1.0, zero_ordinate}};
  std::vector<Weight> expected_weights{1.0};

  EXPECT_THAT(test_scalar_quad_.quadrature_points(),
              ContainerEq(expected_quadrature_points));


  EXPECT_THAT(test_scalar_quad_.quadrature_points_map(),
              ContainerEq(expected_quad_pt_map));

  EXPECT_THAT(test_scalar_quad_.quadrature_weights(),
              ContainerEq(expected_weights));
}


} // namespace