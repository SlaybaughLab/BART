#include "solver/eigenvalue/krylov_schur_eigenvalue_solver.hpp"

#include <deal.II/lac/petsc_full_matrix.h>

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;
using ::testing::ContainerEq;

class SolverEigenvalueKrylovSchurEigenvalueSolverTest : public ::testing::Test {
 public:
  using Matrix = dealii::PETScWrappers::FullMatrix;
  Matrix test_matrix_;
  auto SetUp() -> void override;
};

auto SolverEigenvalueKrylovSchurEigenvalueSolverTest::SetUp() -> void {
  const std::vector<std::vector<double>> matrix_values{ {3, 2, 4}, {2, 0, 2}, {4, 2, 3}};
  test_matrix_.reinit(3,3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      test_matrix_.set(i, j, matrix_values.at(i).at(j));
    }
  }
  test_matrix_.compress(dealii::VectorOperation::insert);
}

TEST_F(SolverEigenvalueKrylovSchurEigenvalueSolverTest, Dummy) {
  std::vector<double> expected_eigenvector{ 2.0/3.0, 1.0/3.0, 2.0/3.0};
  const double expected_spectral_radius{ 8 };
  solver::eigenvalue::KrylovSchurEigenvalueSolver solver;
  const auto [spectral_radius, eigenvector] = solver.SpectralRadius(&test_matrix_);
  EXPECT_NEAR(spectral_radius, expected_spectral_radius, 1e-6);
  ASSERT_EQ(eigenvector.size(), 3);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(std::abs(eigenvector.at(i)), expected_eigenvector.at(i), 1e-6);
  }

}

} // namespace
