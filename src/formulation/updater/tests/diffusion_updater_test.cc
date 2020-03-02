#include "formulation/updater/diffusion_updater.h"

#include "formulation/scalar/tests/diffusion_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationUpdaterDiffusionTest :
    public bart::formulation::updater::test_helpers::UpdaterTests<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(FormulationUpdaterDiffusionTest, bart::testing::AllDimensions);

} // namespace