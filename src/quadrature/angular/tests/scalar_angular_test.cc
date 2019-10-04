#include "quadrature/angular/scalar_angular.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureAngularScalarAngularTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureAngularScalarAngularTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureAngularScalarAngularTest, Order) {
  constexpr int dim = this->dim;
  quadrature::angular::ScalarAngular<dim> test_quad;
  EXPECT_EQ(0, test_quad.order());
}

TYPED_TEST(QuadratureAngularScalarAngularTest, GenerateQuadrature) {
  constexpr int dim = this->dim;
  quadrature::angular::ScalarAngular<dim> test_quad;
  auto quadrature_set = test_quad.GenerateSet();

  ASSERT_EQ(quadrature_set.size(), 1);
  auto quad_point = quadrature_set.at(0);
  auto& [position, weight] = quad_point;

  for (int i = 0; i < dim; ++i)
    EXPECT_EQ(position.get().at(i), 0);

  EXPECT_EQ(weight.get(), 1.0);

}

} // namespace