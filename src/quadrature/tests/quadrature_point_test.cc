#include "quadrature/quadrature_point.h"

#include "quadrature/tests/ordinate_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadraturePointTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  std::shared_ptr<quadrature::OrdinateI<dim>> ordinate_ptr;

  void SetUp() override;
};

template <typename DimensionWrapper>
void QuadraturePointTest<DimensionWrapper>::SetUp() {
  ordinate_ptr = std::make_shared<quadrature::OrdinateMock<dim>>();
}

TYPED_TEST_CASE(QuadraturePointTest, bart::testing::AllDimensions);

TYPED_TEST(QuadraturePointTest, Constructor) {
  constexpr int dim = this->dim;
  const double weight = 1.45;

  quadrature::QuadraturePoint<dim> default_point;
  EXPECT_EQ(default_point.ordinate(), nullptr);
  EXPECT_EQ(default_point.weight(), 0);

  quadrature::QuadraturePoint<dim> new_point(this->ordinate_ptr,
                                             quadrature::Weight(weight));

  EXPECT_EQ(this->ordinate_ptr.use_count(), 2);
  EXPECT_EQ(new_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(new_point.weight(), weight);
}

TYPED_TEST(QuadraturePointTest, Setters) {
  constexpr int dim = this->dim;
  const double weight = 1.45;

  quadrature::QuadraturePoint<dim> test_point;

  auto second_point_ptr = std::make_shared<quadrature::OrdinateMock<dim>>();
  const double second_weight = 2.1;

  test_point.SetOrdinate(this->ordinate_ptr);
  test_point.SetWeight(quadrature::Weight(weight));

  EXPECT_GT(this->ordinate_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(test_point.weight(), weight);

  test_point.SetOrdinate(second_point_ptr);
  test_point.SetWeight(quadrature::Weight(second_weight));

  EXPECT_GT(second_point_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), second_point_ptr);
  EXPECT_EQ(test_point.weight(), second_weight);

  test_point.SetTo(this->ordinate_ptr, quadrature::Weight(weight));
  EXPECT_GT(this->ordinate_ptr.use_count(), 1);
  EXPECT_EQ(test_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(test_point.weight(), weight);
}


} // namespace