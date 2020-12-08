#include "quadrature/calculators/drift_diffusion_integrated_flux.hpp"

#include "quadrature/calculators/tests/drift_diffusion_integrated_flux_mock.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class DriftDiffusionIntegratedFluxTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(DriftDiffusionIntegratedFluxTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionIntegratedFluxTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
