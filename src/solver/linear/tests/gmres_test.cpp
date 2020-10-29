#include "solver/linear/gmres.h"

#include <deal.II/lac/petsc_full_matrix.h>

#include <gtest/gtest.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/petsc_vector.h>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

class SolverGMRESTest : public ::testing::Test {
 protected:
  using FullMatrix = dealii::PETScWrappers::FullMatrix;
  using Vector = dealii::PETScWrappers::MPI::Vector;
  using Preconditioner = dealii::PETScWrappers::PreconditionerBase;
};

TEST_F(SolverGMRESTest, Constructor) {
  bart::solver::linear::GMRES solver;
  EXPECT_EQ(solver.max_iterations(), 100);
  EXPECT_EQ(solver.convergence_tolerance(), 1e-10);

  EXPECT_EQ(solver.solver_control().max_steps(), 100);
  EXPECT_EQ(solver.solver_control().tolerance(), 1e-10);


  bart::solver::linear::GMRES solver_2(210, 1e-6);
  EXPECT_EQ(solver_2.solver_control().max_steps(), 210);
  EXPECT_EQ(solver_2.solver_control().tolerance(), 1e-6);
}

TEST_F(SolverGMRESTest, SolveTestNoPrecon) {

  std::vector<double> b{5,7,8};
  std::vector<double> x{-15, 8, 2};
  std::vector<std::vector<double>> A = {
      {1, 3, -2}, {3, 5, 6}, {2, 4, 3}
  };

  std::vector<unsigned int> indices{0,1,2};
  std::vector<double> zeroes(3,0);

  Vector petsc_b(MPI_COMM_WORLD, 3, 3);
  petsc_b.set(indices, b);
  petsc_b.compress(dealii::VectorOperation::insert);
  Vector petsc_x(MPI_COMM_WORLD, 3, 3);
  petsc_x.set(indices, zeroes);
  petsc_x.compress(dealii::VectorOperation::insert);

  FullMatrix petsc_A(3,3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      petsc_A.set(i, j, A[i][j]);
    }
  }
  petsc_A.compress(dealii::VectorOperation::insert);

  dealii::PETScWrappers::PreconditionNone no_conditioner(petsc_A);

  bart::solver::linear::GMRES solver(100, 1e-6);
  solver.Solve(&petsc_A, &petsc_x, &petsc_b, &no_conditioner);

  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(petsc_x[i], x[i], 1e-6);
  }
}