#include "formulation/updater/saaf_updater.h"

#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationUpdaterSAAFTest : public ::testing::Test {
 public:
};

TYPED_TEST_SUITE(FormulationUpdaterSAAFTest, bart::testing::AllDimensions);

TYPED_TEST(FormulationUpdaterSAAFTest, Constructor) {
  EXPECT_FALSE(true);
}


} // namespace
