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

  quadrature::QuadraturePoint<dim> new_point(this->ordinate_ptr,
                                             quadrature::Weight(weight));

  EXPECT_EQ(this->ordinate_ptr.use_count(), 2);
  EXPECT_EQ(new_point.ordinate(), this->ordinate_ptr);
  EXPECT_EQ(new_point.weight(), weight);
}

} // namespace