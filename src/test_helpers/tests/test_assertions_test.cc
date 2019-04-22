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
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

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
  vector_1 = btest::RandomVector(5, 0, 1.0);
  vector_2 = btest::RandomVector(5, 1.0, 2.0);
  dealii_vector_1.reinit(5);
  dealii_vector_2.reinit(5);

  for (int i = 0; i < 5; ++i) {
    dealii_vector_1[i] = vector_1[i];
    dealii_vector_2[i] = vector_2[i];
  }
}

TEST_F(TestAssertionsTest, GoodComparisonStdVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(vector_1, vector_1));
}

TEST_F(TestAssertionsTest, BadComparisonStdVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(vector_1, vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(vector_2, vector_1));
}

TEST_F(TestAssertionsTest, GoodComparisonDealiiVector) {
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(dealii_vector_1, dealii_vector_1));
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareVector(dealii_vector_2, dealii_vector_2));
}

TEST_F(TestAssertionsTest, BadComparisonDealiiVector) {
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(dealii_vector_1, dealii_vector_2));
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareVector(dealii_vector_2, dealii_vector_1));
}

class TestAssertionsMPIMatricesTests : public::testing::Test {
 protected:
  using Cell = typename dealii::DoFHandler<2>::active_cell_iterator;
  TestAssertionsMPIMatricesTests();
  void SetUp() override;

  dealii::ConstraintMatrix constraint_matrix_;
  dealii::parallel::distributed::Triangulation<2> triangulation_;
  dealii::DoFHandler<2> dof_handler_;
  dealii::FE_Q<2> fe_;
  dealii::IndexSet locally_relevant_dofs;
  dealii::IndexSet locally_owned_dofs_;
  dealii::PETScWrappers::MPI::SparseMatrix matrix_1, matrix_2;
  std::vector<Cell> cells_;
};

TestAssertionsMPIMatricesTests::TestAssertionsMPIMatricesTests()
    : triangulation_(MPI_COMM_WORLD,
                     typename dealii::Triangulation<2>::MeshSmoothing(
                         dealii::Triangulation<2>::smoothing_on_refinement |
                             dealii::Triangulation<2>::smoothing_on_coarsening)),
      dof_handler_(triangulation_),
      fe_(1) {}

void TestAssertionsMPIMatricesTests::SetUp() {
  dealii::GridGenerator::hyper_cube(triangulation_, 0, 1);
  triangulation_.refine_global(2);
  dof_handler_.distribute_dofs(fe_);
  dealii::DoFTools::extract_locally_relevant_dofs(dof_handler_,
                                                  locally_relevant_dofs);
  locally_owned_dofs_ = dof_handler_.locally_owned_dofs();

  for (auto cell = dof_handler_.begin_active(); cell != dof_handler_.end(); ++ cell) {
    if (cell->is_locally_owned())
      cells_.push_back(cell);
  }

  constraint_matrix_.clear();
  constraint_matrix_.reinit(locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler_,
                                                  constraint_matrix_);
  constraint_matrix_.close();

  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler_, dsp,
                                          constraint_matrix_, false);
  dealii::SparsityTools::distribute_sparsity_pattern(
      dsp,
      dof_handler_.n_locally_owned_dofs_per_processor(),
      MPI_COMM_WORLD, locally_relevant_dofs);
  constraint_matrix_.condense(dsp);

  matrix_1.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp, MPI_COMM_WORLD);
  matrix_2.reinit(locally_owned_dofs_, locally_owned_dofs_, dsp, MPI_COMM_WORLD);

  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      for (auto index_j : local_dof_indices) {
        matrix_1.add(index_i, index_j, 1);
        matrix_2.add(index_i, index_j, 2);
      }
    }
    matrix_1.compress(dealii::VectorOperation::add);
    matrix_2.compress(dealii::VectorOperation::add);
  }
}

TEST_F(TestAssertionsMPIMatricesTests, CompareMPIMatrices) {
  EXPECT_EQ(AssertionSuccess(),
            bart::testing::CompareMPIMatrices(matrix_1, matrix_1));
  EXPECT_EQ(AssertionFailure(),
            bart::testing::CompareMPIMatrices(matrix_1, matrix_2));
}

} // namespace

