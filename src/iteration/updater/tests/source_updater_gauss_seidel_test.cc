#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>

#include "test_helpers/gmock_wrapper.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "formulation/cfem_stamper_i.h"

namespace  {

using namespace bart;

class IterationSourceUpdaterGaussSeidelTest : public ::testing::Test {
 protected:
  using CFEMSourceUpdater =
      iteration::updater::SourceUpdaterGaussSeidel<formulation::CFEMStamperI>;
  std::unique_ptr<formulation::CFEMStamperI> mock_stamper_ptr_;

  void SetUp() override;
};

void IterationSourceUpdaterGaussSeidelTest::SetUp() {
  mock_stamper_ptr_ = std::make_unique<formulation::CFEM_StamperMock>();
}



TEST_F(IterationSourceUpdaterGaussSeidelTest, Constructor) {
  CFEMSourceUpdater test_updater(std::move(mock_stamper_ptr_));
  EXPECT_EQ(mock_stamper_ptr_, nullptr);
}

} // namespace