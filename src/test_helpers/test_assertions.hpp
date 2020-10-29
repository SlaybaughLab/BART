#ifndef BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_
#define BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart::test_helpers {

using FullMatrix = dealii::FullMatrix<double>;
using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
using MPIVector = dealii::PETScWrappers::MPI::Vector;

using ::testing::AssertionResult;

[[nodiscard]] auto AreEqual(const dealii::Vector<double>& expected, const dealii::Vector<double>& result,
                            const double tol = 1e-6) -> AssertionResult;

[[nodiscard]] auto AreEqual(const std::vector<double> expected, const std::vector<double> result,
                            const double tol = 1e-6) -> AssertionResult;

[[nodiscard]] auto AreEqual(const FullMatrix& expected, const FullMatrix& result,
                            const double tol = 1e-6) -> AssertionResult;

[[nodiscard]] auto AreEqual(const MPISparseMatrix& expected, const MPISparseMatrix& result) -> AssertionResult;

[[nodiscard]] auto CompareMPIVectors(const MPIVector& expected, const MPIVector& result) -> AssertionResult;

} // namespace bart::test_helpers

#endif // BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_