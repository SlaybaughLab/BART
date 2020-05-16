#include "test_helpers/test_assertions.h"

#include <vector>

#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/petsc_sparse_matrix.h>

#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::AssertionResult;
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;

class TestAssertionsTest : public ::testing::Test {
 protected:
  std::vector<double> vector_1;
  std::vector<double> vector_2;
  dealii::Vector<double> dealii_vector_1;
  dealii::Vector<double> dealii_vector_2;

  void SetUp() override;
};

void TestAssertionsTest::SetUp() {
  vector_1 = test_helpers::RandomVector(5, 0, 1.0);
  vector_2 = test_helpers::RandomVector(5, 1.0, 2.0);
  dealii_vector_1.reinit(5);
  dealii_vector_2.reinit(5);

  for (int i = 0; i < 5; ++i) {
    dealii_vector_1[i] = vector_1[i];
    dealii_vector_2[i] = vector_2[i];
  }
}

TEST_F(TestAssertionsTest, GoodComparisonStdVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareVector(vector_1, vector_1));
}

TEST_F(TestAssertionsTest, BadComparisonStdVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareVector(vector_1, vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareVector(vector_2, vector_1));
}

TEST_F(TestAssertionsTest, GoodComparisonDealiiVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareVector(dealii_vector_1, dealii_vector_1));
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareVector(dealii_vector_2, dealii_vector_2));
}

TEST_F(TestAssertionsTest, BadComparisonDealiiVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareVector(dealii_vector_1, dealii_vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareVector(dealii_vector_2, dealii_vector_1));
}

class TestAssertionsMPIMatricesTests : public ::testing::Test,
                                       public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
};

void TestAssertionsMPIMatricesTests::SetUp() {
  SetUpDealii();
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      for (auto index_j : local_dof_indices) {
        matrix_1.add(index_i, index_j, 1);
        matrix_2.add(index_i, index_j, 2);
        matrix_3.add(index_i, index_j, 1);
      }
    }
  }

  matrix_1.compress(dealii::VectorOperation::add);
  matrix_2.compress(dealii::VectorOperation::add);
  matrix_3.compress(dealii::VectorOperation::add);
}

TEST_F(TestAssertionsMPIMatricesTests, CompareMPIMatrices) {
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareMPIMatrices(matrix_2, matrix_2));
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareMPIMatrices(matrix_1, matrix_3));
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareMPIMatrices(matrix_1, matrix_2));

  int random_cell = test_helpers::RandomDouble(0, cells_.size());
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);
  cells_[random_cell]->get_dof_indices(local_dof_indices);

  for (auto index_i : local_dof_indices) {
    for (auto index_j : local_dof_indices) {
      matrix_3.add(index_i, index_j, 1);
    }
  }

  matrix_3.compress(dealii::VectorOperation::add);

  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareMPIMatrices(matrix_1, matrix_3));
}

class TestAssertionsMPIVectorTests : public ::testing::Test,
                                     public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
};

void TestAssertionsMPIVectorTests::SetUp() {
  SetUpDealii();
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      vector_1(index_i) += 1;
      vector_2(index_i) += 2;
      vector_3(index_i) += 1;
    }
  }

  vector_1.compress(dealii::VectorOperation::add);
  vector_2.compress(dealii::VectorOperation::add);
  vector_3.compress(dealii::VectorOperation::add);
}

TEST_F(TestAssertionsMPIVectorTests, CompareMPIVectors) {
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareMPIVectors(vector_1, vector_1));
  EXPECT_EQ(AssertionSuccess(),
            bart::test_helpers::CompareMPIVectors(vector_1, vector_3));
  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareMPIVectors(vector_1, vector_2));

  int random_cell = test_helpers::RandomDouble(0, cells_.size());
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);
  cells_[random_cell]->get_dof_indices(local_dof_indices);

  for (auto index_i : local_dof_indices) {
    vector_3(index_i) += 1;
  }

  EXPECT_EQ(AssertionFailure(),
            bart::test_helpers::CompareMPIVectors(vector_1, vector_3));
}

} // namespace

