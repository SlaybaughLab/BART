#include "iteration/updater/fixed_updater.h"

#include <memory>

#include <deal.II/lac/full_matrix.h>

#include "data/system.h"
#include "data/system/tests/bilinear_term_mock.h"
#include "formulation/tests/cfem_stamper_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"

namespace {

using namespace bart;

using ::testing::An, ::testing::NiceMock, ::testing::Return;

/*
 * ===== BASIC TESTS ===========================================================
 * Tests operation of the constructor, provides mock stamper and observing ptr
 * to stamper.
 *
 */

class IterationFixedUpdaterBasicTest : public ::testing::Test {
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

void IterationFixedUpdaterBasicTest::SetUp() {
  mock_stamper_ptr_ = std::make_unique<NiceMock<MockStamperType>>();
  test_updater_ptr_ = std::make_unique<FixedUpdaterType>(std::move(mock_stamper_ptr_));
  stamper_obs_ptr_ = test_updater_ptr_->GetStamperPtr();
}

TEST_F(IterationFixedUpdaterBasicTest, Constructor) {
  // Constructor should have stored passed stamper
  EXPECT_FALSE(stamper_obs_ptr_ == nullptr);

  // Constructor should throw an error if Stamper object is invalid
  std::unique_ptr<MockStamperType> empty_stamper_ptr_ ;
  EXPECT_ANY_THROW({
    FixedUpdaterType throwing_updater(std::move(empty_stamper_ptr_));
  });
}

/* ===== DOMAIN TESTS ==========================================================
 * Tests updating of system matrices. DealiiTestDomain is used to provide
 * matrices with sparsity patterns that can be stamped properly. A test system
 * struct will be instantiated (populated with mock classes) and passed to
 * the updater to verify it mediates actions between the system and the
 * stamper correctly.
 */


class IterationFixedUpdaterDomainTest : public IterationFixedUpdaterBasicTest,
                                        public bart::testing::DealiiTestDomain<2> {
 protected:

  data::System test_system_;
  std::unique_ptr<data::system::BilinearTermMock> mock_lhs_ptr_;
  std::shared_ptr<data::system::MPISparseMatrix> matrix_ptr_;

  void SetUp() override;
};

void IterationFixedUpdaterDomainTest::SetUp() {
  IterationFixedUpdaterBasicTest::SetUp();
  SetUpDealii();
  // Stamp matrices with specified values
  StampMatrix(matrix_1, 1);

  // This mock left hand side will be stored in our system and will return
  // our class matrix_ptr_ by default.
  mock_lhs_ptr_ = std::make_unique<NiceMock<data::system::BilinearTermMock>>();

  // Setup matrix_ptr to make it identical to the DealiiTestDomain matrices,
  // then stamp with the value 2
  matrix_ptr_ = std::make_shared<data::system::MPISparseMatrix>();
  matrix_ptr_->reinit(matrix_1);
  StampMatrix(*matrix_ptr_, 2);

  ON_CALL(*mock_lhs_ptr_, GetFixedTermPtr(An<data::system::Index>()))
      .WillByDefault(Return(matrix_ptr_));
  ON_CALL(*mock_lhs_ptr_, GetFixedTermPtr(An<data::system::GroupNumber>()))
      .WillByDefault(Return(matrix_ptr_));
}

TEST_F(IterationFixedUpdaterDomainTest, Dummy) {
  EXPECT_FALSE(bart::testing::CompareMPIMatrices(matrix_1, matrix_2));
}



} // namespace
