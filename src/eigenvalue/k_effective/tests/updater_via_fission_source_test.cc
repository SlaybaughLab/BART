#include <calculator/cell/tests/integrated_fission_source_mock.h>
#include "eigenvalue/k_effective/updater_via_fission_source.h"

#include "calculator/cell/tests/total_aggregated_fission_source_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class EigKEffUpdaterViaFissionSourceTest : public ::testing::Test {
 protected:

  using TotalAggregatedFissionSourceType = calculator::cell::TotalAggregatedFissionSourceMock;

  // Supporting objects
  std::unique_ptr<TotalAggregatedFissionSourceType> fission_source_mock_ptr_;

  void SetUp() override;
};

void EigKEffUpdaterViaFissionSourceTest::SetUp() {
  fission_source_mock_ptr_ =
      std::make_unique<TotalAggregatedFissionSourceType>();
}

TEST_F(EigKEffUpdaterViaFissionSourceTest, Constructor) {
  eigenvalue::k_effective::UpdaterViaFissionSource test_k_eff_updater;

  EXPECT_EQ(test_k_eff_updater.k_effective(), 1);
  EXPECT_EQ(test_k_eff_updater.previous_fission_source(), 0);
  EXPECT_EQ(test_k_eff_updater.current_fission_source(), 0);
}

} // namespace