#ifndef BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_
#define BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_

#include <deal.II/lac/vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace testing {

using ::testing::AssertionResult;
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;

AssertionResult CompareVector(const dealii::Vector<double>& expected,
                              const dealii::Vector<double>& result,
                              const double tol = 1e-6);

AssertionResult CompareVector(const std::vector<double> expected,
                              const std::vector<double> result,
                              const double tol = 1e-6);

AssertionResult CompareMPIMatrices(
    const dealii::PETScWrappers::MPI::SparseMatrix& expected,
    const dealii::PETScWrappers::MPI::SparseMatrix& result);

} // namespace testing

} // namespace bart

#endif // BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_