#include "test_helpers/gmock_wrapper.h"

#include "calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DriftDiffusionVectorCalculatorTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(DriftDiffusionVectorCalculatorTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionVectorCalculatorTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
