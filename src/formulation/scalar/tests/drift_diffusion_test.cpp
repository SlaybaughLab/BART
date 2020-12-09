#include "formulation/scalar/drift_diffusion.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DriftDiffusionFormulationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(DriftDiffusionFormulationTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionFormulationTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
