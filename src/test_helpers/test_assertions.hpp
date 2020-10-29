#ifndef BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_
#define BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_

#include <deal.II/lac/vector.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart::test_helpers {

using MPISparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
using MPIVector = dealii::PETScWrappers::MPI::Vector;

using ::testing::AssertionResult;

[[nodiscard]] auto CompareVector(const dealii::Vector<double>& expected,
                                            const dealii::Vector<double>& result,
                                            const double tol = 1e-6) -> AssertionResult;

[[nodiscard]] auto CompareVector(const std::vector<double> expected, const std::vector<double> result,
                                            const double tol = 1e-6) -> AssertionResult;

AssertionResult CompareMPIMatrices(
    const dealii::PETScWrappers::MPI::SparseMatrix& expected,
    const dealii::PETScWrappers::MPI::SparseMatrix& result);

[[nodiscard]] auto CompareMPIMatrices(const MPISparseMatrix& expected,
                                      const MPISparseMatrix& result) -> AssertionResult;

[[nodiscard]] auto CompareMPIVectors(const MPIVector& expected, const MPIVector& result) -> AssertionResult;

} // namespace bart::test_helpers

#endif // BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_HPP_