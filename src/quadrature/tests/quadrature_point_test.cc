#include "quadrature/quadrature_point.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadraturePointTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadraturePointTest, bart::testing::AllDimensions);

TYPED_TEST(QuadraturePointTest, Constructor) {
  EXPECT_TRUE(true);
}

} // namespace