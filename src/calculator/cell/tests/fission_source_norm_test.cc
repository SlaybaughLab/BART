#include "calculator/cell/fission_source_norm.h"

#include "test_helpers/gmock_wrapper.h"


namespace  {

using namespace bart;

template <typename DimensionWrapper>
class CalcCellFissionSourceNormTest :public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(CalcCellFissionSourceNormTest, bart::testing::AllDimensions);

TYPED_TEST(CalcCellFissionSourceNormTest, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace