#include "iteration/updater/source_updater_gauss_seidel.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IterationSourceUpdaterGaussSeidelTest : public ::testing::Test {
 protected:
  iteration::updater::SourceUpdaterGaussSeidel test_updater;
};

TEST_F(IterationSourceUpdaterGaussSeidelTest, Dummy) {
  EXPECT_TRUE(true);
}

} // namespace