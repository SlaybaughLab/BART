#include "quadrature/angular/angular_quadrature_scalar.h"

#include "quadrature/angular/angular_quadrature_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

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
  auto dim = this->dim;
  EXPECT_EQ(test_scalar_quad_.total_quadrature_points(), 1);
}


} // namespace