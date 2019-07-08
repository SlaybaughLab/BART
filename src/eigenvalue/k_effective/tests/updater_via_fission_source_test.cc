#include "eigenvalue/k_effective/updater_via_fission_source.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

 class EigKEffUpdaterViaFissionSourceTest : public ::testing::Test {
 protected:
};

TEST_F(EigKEffUpdaterViaFissionSourceTest, Constructor) {
  eigenvalue::k_effective::UpdaterViaFissionSource test_k_eff_updater;

  EXPECT_EQ(test_k_eff_updater.k_effective(), 1);
  EXPECT_EQ(test_k_eff_updater.previous_fission_source(), 0);
  EXPECT_EQ(test_k_eff_updater.current_fission_source(), 0);
}

} // namespace