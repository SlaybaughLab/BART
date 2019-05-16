#include "iteration/initializer/set_fixed_terms.h"

#include "data/system.h"
#include "data/system/tests/bilinear_term_mock.h"
#include "iteration/updater/tests/fixed_updater_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::NiceMock;

class IterationInitializerSetFixedTermsTest : public ::testing::Test {
 protected:
  using InitializerType = iteration::initializer::SetFixedTerms;
  using MockFixedUpdaterType = iteration::updater::FixedUpdaterMock;
  using MockBilinearTermType = data::system::BilinearTermMock;

  // Initializer to be tested
  std::unique_ptr<InitializerType> test_initializer_;

  // Supporting objects
  data::System test_system_;
  data::system::MPISparseMatrix mpi_matrix_;

  // Pointers for observing mocks owned by other objects
  MockBilinearTermType* bilinear_term_obs_ptr_ = nullptr;
  MockFixedUpdaterType* updater_obs_ptr_ = nullptr;

  void SetUp() override;
};

void IterationInitializerSetFixedTermsTest::SetUp() {
  // Set up testing object
  auto mock_fixed_updater_ptr =
      std::make_unique<NiceMock<MockFixedUpdaterType>>();
  test_initializer_ =
      std::make_unique<InitializerType>(std::move(mock_fixed_updater_ptr));

  // Set up supporting objects
  auto mock_bilinear_term = std::make_unique<MockBilinearTermType>();
  test_system_.left_hand_side_ptr_ = std::move(mock_bilinear_term);

  // Set up observing pointers
  updater_obs_ptr_ =
      dynamic_cast<MockFixedUpdaterType*>(test_initializer_->GetUpdater());
  bilinear_term_obs_ptr_ =
      dynamic_cast<MockBilinearTermType*>(test_system_.left_hand_side_ptr_.get());
}

TEST_F(IterationInitializerSetFixedTermsTest, Constructor) {
  EXPECT_TRUE(updater_obs_ptr_ != nullptr);
}

} // namespace