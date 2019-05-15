#include "iteration/updater/fixed_updater.h"

#include <memory>

#include "formulation/tests/cfem_stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::NiceMock;

class IterationFixedUpdaterTest : public ::testing::Test {
 protected:
  // Stamper interface type
  using StamperType = formulation::CFEMStamperI;
  // Mock stamper object type (derived from StamperType)
  using MockStamperType = formulation::CFEM_StamperMock;
  // Fixed updater type (based on stamper interface)
  using FixedUpdaterType = iteration::updater::FixedUpdater<StamperType>;

  std::unique_ptr<MockStamperType> mock_stamper_ptr_;
  std::unique_ptr<FixedUpdaterType> test_updater_ptr_;
  StamperType* stamper_obs_ptr_ = nullptr;

  void SetUp() override;
};

void IterationFixedUpdaterTest::SetUp() {
  mock_stamper_ptr_ = std::make_unique<NiceMock<MockStamperType>>();
  test_updater_ptr_ = std::make_unique<FixedUpdaterType>(std::move(mock_stamper_ptr_));
  stamper_obs_ptr_ = test_updater_ptr_->GetStamperPtr();
}

TEST_F(IterationFixedUpdaterTest, Constructor) {
  // Constructor should have stored passed stamper
  EXPECT_FALSE(stamper_obs_ptr_ == nullptr);

  // Constructor should throw an error if Stamper object is invalid
  std::unique_ptr<MockStamperType> empty_stamper_ptr_ ;
  EXPECT_ANY_THROW({
    FixedUpdaterType throwing_updater(std::move(empty_stamper_ptr_));
  });

}

} // namespace
