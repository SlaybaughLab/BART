#include "test_helpers/test_assertions.hpp"

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
#include <deal.II/base/tensor.h>

#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace test_helpers = bart::test_helpers;

using ::testing::AssertionResult, ::testing::AssertionFailure, ::testing::AssertionSuccess;

class TestAssertionsTest : public ::testing::Test {
 protected:
  std::vector<double> vector_1, vector_2;
  dealii::Vector<double> dealii_vector_1, dealii_vector_2;
  void SetUp() override;
};

void TestAssertionsTest::SetUp() {
  const auto vector_size = test_helpers::RandomInt(5, 10);
  vector_1 = test_helpers::RandomVector(vector_size, 0, 1.0);
  vector_2 = test_helpers::RandomVector(vector_size, 1.0, 2.0);
  dealii_vector_1.reinit(vector_size);
  dealii_vector_2.reinit(vector_size);

  for (int i = 0; i < vector_size; ++i) {
    dealii_vector_1[i] = vector_1[i];
    dealii_vector_2[i] = vector_2[i];
  }
}

TEST_F(TestAssertionsTest, GoodComparisonStdVector) {
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(vector_1, vector_1));
}

TEST_F(TestAssertionsTest, BadComparisonStdVector) {
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(vector_1, vector_2));
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(vector_2, vector_1));
}

TEST_F(TestAssertionsTest, GoodComparisonDealiiVector) {
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(dealii_vector_1, dealii_vector_1));
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(dealii_vector_2, dealii_vector_2));
}

TEST_F(TestAssertionsTest, BadComparisonDealiiVector) {
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(dealii_vector_1, dealii_vector_2));
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(dealii_vector_2, dealii_vector_1));
}

class TestAssertionsMatrixTests : public ::testing::Test {
 public:
  dealii::FullMatrix<double> matrix_1, matrix_2, matrix_bad_columns, matrix_bad_rows;
  void SetUp() override;
};

void TestAssertionsMatrixTests::SetUp() {
  const auto matrix_rows{ test_helpers::RandomInt(5, 10) }, matrix_cols{ matrix_rows + 1 };
  matrix_1.reinit(matrix_rows, matrix_cols);
  matrix_2.reinit(matrix_rows, matrix_cols);
  matrix_bad_columns.reinit(matrix_rows, matrix_rows);
  matrix_bad_rows.reinit(matrix_cols, matrix_cols);

  for (int i = 0; i < matrix_rows; ++i) {
    for (int j = 0; j < matrix_cols; ++j) {
      matrix_1.set(i, j, test_helpers::RandomDouble(-100, 100));
      matrix_2.set(i, j, test_helpers::RandomDouble(-100, 100));
    }
  }

  for (int i = 0; i < matrix_rows; ++i) {
    for (int j = 0; j < matrix_rows; ++j) {
      matrix_bad_columns.set(i, j, test_helpers::RandomDouble(-100, 100));
    }
  }

  for (int i = 0; i < matrix_cols; ++i) {
    for (int j = 0; j < matrix_cols; ++j) {
      matrix_bad_rows.set(i, j, test_helpers::RandomDouble(-100, 100));
    }
  }
}

TEST_F(TestAssertionsMatrixTests, GoodComparison) {
  EXPECT_TRUE(test_helpers::AreEqual(matrix_1, matrix_1));
  EXPECT_TRUE(test_helpers::AreEqual(matrix_2, matrix_2));
  EXPECT_TRUE(test_helpers::AreEqual(matrix_bad_columns, matrix_bad_columns));
  EXPECT_TRUE(test_helpers::AreEqual(matrix_bad_rows, matrix_bad_rows));
}

TEST_F(TestAssertionsMatrixTests, BadComparison) {
  EXPECT_FALSE(test_helpers::AreEqual(matrix_1, matrix_2));
  EXPECT_FALSE(test_helpers::AreEqual(matrix_2, matrix_1));
}

TEST_F(TestAssertionsMatrixTests, BadSizeComparison) {
  EXPECT_FALSE(test_helpers::AreEqual(matrix_1, matrix_bad_columns));
  EXPECT_FALSE(test_helpers::AreEqual(matrix_2, matrix_bad_columns));
  EXPECT_FALSE(test_helpers::AreEqual(matrix_1, matrix_bad_rows));
  EXPECT_FALSE(test_helpers::AreEqual(matrix_2, matrix_bad_rows));
  EXPECT_FALSE(test_helpers::AreEqual(matrix_bad_rows, matrix_bad_columns));
}

TEST_F(TestAssertionsMatrixTests, GoodComparisonWithinTolerance) {
  auto matrix_3 = matrix_1;
  for (auto entry : matrix_3) {
    entry += test_helpers::RandomDouble(1e-6, 1e-5);
  }
  EXPECT_FALSE(test_helpers::AreEqual(matrix_1, matrix_3));
  EXPECT_TRUE(test_helpers::AreEqual(matrix_1, matrix_3, 1e-4));
}

template <typename DimensionWrapper>
class TestAssertionsTensorsAreEqual : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  dealii::Tensor<1, dim> tensor_1_, tensor_2_;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto TestAssertionsTensorsAreEqual<DimensionWrapper>::SetUp() -> void {
  for (int i = 0; i < dim; ++i) {
    tensor_1_[i] = test_helpers::RandomDouble(-100, 100);
    tensor_2_[i] = test_helpers::RandomDouble(-100, 100);
  }
}

TYPED_TEST_SUITE(TestAssertionsTensorsAreEqual, bart::testing::AllDimensions);

TYPED_TEST(TestAssertionsTensorsAreEqual, GoodComparison) {
  EXPECT_TRUE(test_helpers::AreEqual(this->tensor_1_, this->tensor_1_));
  EXPECT_TRUE(test_helpers::AreEqual(this->tensor_2_, this->tensor_2_));
}

TYPED_TEST(TestAssertionsTensorsAreEqual, BadComparison) {
  EXPECT_FALSE(test_helpers::AreEqual(this->tensor_1_, this->tensor_2_));
  EXPECT_FALSE(test_helpers::AreEqual(this->tensor_2_, this->tensor_1_));
}

TYPED_TEST(TestAssertionsTensorsAreEqual, Tolerance) {
  auto tensor_3 = this->tensor_1_;
  tensor_3 *= (1 + 1e-8);
  EXPECT_TRUE(test_helpers::AreEqual(this->tensor_1_, tensor_3));
  tensor_3 *= (1 + 1e-5);
  EXPECT_FALSE(test_helpers::AreEqual(this->tensor_1_, tensor_3));
}


class TestAssertionsMPIMatricesTests : public ::testing::Test, public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
};

void TestAssertionsMPIMatricesTests::SetUp() {
  SetUpDealii();
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (const auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (const auto index_i : local_dof_indices) {
      for (const auto index_j : local_dof_indices) {
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
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(matrix_2, matrix_2));
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(matrix_1, matrix_3));
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(matrix_1, matrix_2));

  const int random_cell = test_helpers::RandomInt(0, cells_.size());
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);
  cells_[random_cell]->get_dof_indices(local_dof_indices);

  for (const auto index_i : local_dof_indices) {
    for (const auto index_j : local_dof_indices) {
      matrix_3.add(index_i, index_j, 1);
    }
  }
  matrix_3.compress(dealii::VectorOperation::add);
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(matrix_1, matrix_3));
}

class TestAssertionsMPIVectorTests : public ::testing::Test, public bart::testing::DealiiTestDomain<2> {
 protected:
  void SetUp() override;
};

void TestAssertionsMPIVectorTests::SetUp() {
  SetUpDealii();
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (const auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (const auto index_i : local_dof_indices) {
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
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(vector_1, vector_1));
  EXPECT_EQ(AssertionSuccess(), bart::test_helpers::AreEqual(vector_1, vector_3));
  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(vector_1, vector_2));

  const int random_cell = test_helpers::RandomInt(0, cells_.size());
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);
  cells_[random_cell]->get_dof_indices(local_dof_indices);

  for (auto index_i : local_dof_indices) {
    vector_3(index_i) += 1;
  }

  EXPECT_EQ(AssertionFailure(), bart::test_helpers::AreEqual(vector_1, vector_3));
}

} // namespace

