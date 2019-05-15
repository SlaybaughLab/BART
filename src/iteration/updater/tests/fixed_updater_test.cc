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

using ::testing::NiceMock;

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
 * matrices with sparsity patterns that can be stamped properly.
 */


class IterationFixedUpdaterDomainTest : public IterationFixedUpdaterBasicTest,
                                        public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
  void StampMatrix(data::system::MPISparseMatrix& to_stamp, double value);
};

void IterationFixedUpdaterDomainTest::SetUp() {
  IterationFixedUpdaterBasicTest::SetUp();
  SetUpDealii();
  StampMatrix(matrix_1, 1);
  StampMatrix(matrix_2, 2);
}
// Stamps a matrix with a given value or 1
void IterationFixedUpdaterDomainTest::StampMatrix(
    data::system::MPISparseMatrix &to_stamp,
    const double value = 1) {
  dealii::FullMatrix<double> cell_matrix(fe_.dofs_per_cell, fe_.dofs_per_cell);

  for (unsigned int i = 0; i < cell_matrix.m(); ++i) {
    for (unsigned int j = 0; j < cell_matrix.n(); ++j) {
      cell_matrix(i,j) = value;
    }
  }

  for (const auto& cell : cells_) {
    std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
    to_stamp.add(local_dof_indices, cell_matrix);
  }

  to_stamp.compress(dealii::VectorOperation::add);
}

TEST_F(IterationFixedUpdaterDomainTest, Dummy) {
  EXPECT_FALSE(bart::testing::CompareMPIMatrices(matrix_1, matrix_2));
}



} // namespace
