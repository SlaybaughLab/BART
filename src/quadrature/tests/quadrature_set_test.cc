#include "quadrature/quadrature_set.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureSetTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureSetTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureSetTest, Constructor) {
  constexpr int dim = this->dim;
  EXPECT_NO_THROW(quadrature::QuadratureSet<dim> test_set);
}

} // namespace