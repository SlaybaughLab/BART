#include "quadrature/angular/angular_quadrature_scalar.h"

#include "quadrature/angular/angular_quadrature_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class AngularQuadratureScalarTests : public ::testing::Test {
 protected:
  // Tested object
  quadrature::angular::AngularQuadratureScalar test_scalar_quad_;
};

TEST_F(AngularQuadratureScalarTests, Dummy) {
  EXPECT_TRUE(true);
}


} // namespace