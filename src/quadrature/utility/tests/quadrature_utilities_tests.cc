#include "quadrature/utility/quadrature_utilities.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class QuadratureUtilityTests : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureUtilityTests, bart::testing::AllDimensions);

TYPED_TEST(QuadratureUtilityTests, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace