#include "test_helpers/test_assertions.hpp"

#include <deal.II/base/mpi.h>

namespace bart {

namespace test_helpers {

namespace {
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;

using DealiiVector = dealii::Vector<double>;
using FullMatrix = dealii::FullMatrix<double>;
}

template <>
auto AreEqual<dealii::Vector<double>>(const DealiiVector& expected, const DealiiVector& result,
                                      const double tol) -> AssertionResult {
  unsigned int size = expected.size();

  if (result.size() != size)
    return AssertionFailure() << "Result has wrong number of entries: "
                              << result.size() << ", expected" << size;

  for (unsigned int i = 0; i < size; ++i) {
    if (abs(result[i] - expected[i]) > tol) {
      return AssertionFailure() << "Entry (" << i <<
                                ") has value: " << result[i] <<
                                ", expected: " << expected[i];
    }
  }
  return AssertionSuccess();
}

template <>
auto AreEqual<std::vector<double>>(const std::vector<double>& expected, const std::vector<double>& result,
                                   const double tol) -> AssertionResult{
  DealiiVector expected_vec(expected.begin(), expected.end());
  DealiiVector result_vec(result.begin(), result.end());
  return AreEqual(expected_vec, result_vec, tol);
}

template <>
auto AreEqual<FullMatrix>(const FullMatrix &expected, const FullMatrix &result,
                          const double tol) -> AssertionResult {
  using size_type = FullMatrix::size_type;
  const size_type rows{expected.m()}, cols{expected.n()};

  if (const size_type result_rows = result.m(); rows != result_rows) {
    return AssertionFailure() << "Expected matrix has n_rows = " << rows << ", while result matrix has n_rows = "
                              << result_rows;
  } else if (const size_type result_cols = result.n(); cols != result_cols) {
    return AssertionFailure() << "Expected matrix has n_cols = " << cols << ", while result matrix has n_cols = "
                              << result_cols;
  }

  for (size_type i = 0; i < rows; ++i) {
    for (size_type j = 0; j < cols; ++j) {
      if (abs(expected(i, j) - result(i,j)) > tol) {
        return AssertionFailure() << "Expected matrix has value " << expected(i,j) << " at (" << i << ", " << j
                                  << ") while result matrix has " << result(i, j);
      }
    }
  }


  return AssertionSuccess();
}

AssertionResult AreEqual(
    const dealii::PETScWrappers::MPI::SparseMatrix& expected,
    const dealii::PETScWrappers::MPI::SparseMatrix& result) {

  auto [first_local_row, last_local_row] = expected.local_range();
  unsigned int n_columns = expected.n();
  bool has_failed = false;

  int this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  dealii::PETScWrappers::MPI::Vector results(MPI_COMM_WORLD,
                                             n_processes, 1);

  for (unsigned int i = first_local_row; i < last_local_row; ++i) {
    for (unsigned int j = 0; j < n_columns; ++j) {
      if (result(i, j) != expected(i, j)) {
        results(this_process) += 1;
        has_failed = true;
      }
      if (has_failed)
        break;
    }
    if (has_failed)
      break;
  }

  results(this_process) += 0;
  if (results.l1_norm() > 0) {
    return AssertionFailure();
  } else {
    return AssertionSuccess();
  }
}

AssertionResult AreEqual(
    const dealii::PETScWrappers::MPI::Vector& expected,
    const dealii::PETScWrappers::MPI::Vector& result) {

  auto [first_local_row, last_local_row] = expected.local_range();
  bool has_failed = false;

  int this_process = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  int n_processes = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  dealii::PETScWrappers::MPI::Vector results(MPI_COMM_WORLD,
                                             n_processes, 1);

  for (unsigned int i = first_local_row; i < last_local_row; ++i) {
    if (result(i) != expected(i)) {
      results(this_process) += 1;
      has_failed = true;
    }
    if (has_failed)
      break;
  }

  results(this_process) += 0;
  if (results.l1_norm() > 0) {
    return AssertionFailure();
  } else {
    return AssertionSuccess();
  }
}

} // namespace test_helpers

} // namespace bart