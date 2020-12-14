#include "formulation/updater/drift_diffusion_updater.hpp"

#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
template <int dim>
using UpdaterTest = bart::formulation::updater::test_helpers::UpdaterTests<dim>;


template <typename DimensionWrapper>
class FormulationUpdaterDriftDiffusionTest : public UpdaterTest<DimensionWrapper::value> {
 public:
  static constexpr int dim{ DimensionWrapper::value };
};

TYPED_TEST_SUITE(FormulationUpdaterDriftDiffusionTest, bart::testing::AllDimensions);

TYPED_TEST(FormulationUpdaterDriftDiffusionTest, Dummy) {
  EXPECT_TRUE(false);
}

} // namespace
