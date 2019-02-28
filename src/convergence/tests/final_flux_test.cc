#include "convergence/final_flux.h"

#include <memory>
#include <optional>

#include "convergence/flux/tests/multi_checker_mock.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

using ::testing::_;

class ConvergenceFinalFluxTest : public ::testing::Test {
 protected:
};

TEST_F(ConvergenceFinalFluxTest, Dummy) {
  ASSERT_TRUE(true);
}
