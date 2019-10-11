#include "quadrature/quadrature_set.h"

#include <memory>

#include "quadrature/tests/quadrature_point_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureSetTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  quadrature::QuadratureSet<dim> test_set_;
};

TYPED_TEST_CASE(QuadratureSetTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureSetTest, Constructor) {
  constexpr int dim = this->dim;
  EXPECT_NO_THROW(quadrature::QuadratureSet<dim> test_set);
}

TYPED_TEST(QuadratureSetTest, AddPoint) {
  constexpr int dim = this->dim;
  auto quadrature_point =
      std::make_shared<quadrature::QuadraturePointMock<dim>>();
  auto second_quadrature_point =
      std::make_shared<quadrature::QuadraturePointMock<dim>>();

  EXPECT_EQ(this->test_set_.size(), 0);

  EXPECT_TRUE(this->test_set_.AddPoint(quadrature_point));
  EXPECT_EQ(this->test_set_.size(), 1);

  EXPECT_TRUE(this->test_set_.AddPoint(second_quadrature_point));
  EXPECT_EQ(this->test_set_.size(), 2);

  EXPECT_FALSE(this->test_set_.AddPoint(quadrature_point));
  EXPECT_EQ(this->test_set_.size(), 2);
}

} // namespace