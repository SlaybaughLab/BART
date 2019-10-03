#include "quadrature/ordinate.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

template <typename DimensionWrapper>
class QuadratureOrdinateTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};


TYPED_TEST_CASE(QuadratureOrdinateTest, bart::testing::AllDimensions);

TYPED_TEST(QuadratureOrdinateTest, Construction) {

}


} // namespace