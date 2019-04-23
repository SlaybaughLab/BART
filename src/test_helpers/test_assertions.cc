#include "test_helpers/test_assertions.h"

#include <deal.II/base/mpi.h>

namespace bart {

namespace testing {

AssertionResult CompareVector(const dealii::Vector<double>& expected,
                                     const dealii::Vector<double>& result,
                                     const double tol) {
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

AssertionResult CompareVector(const std::vector<double> expected,
                                     const std::vector<double> result,
                                     const double tol) {
  dealii::Vector<double> expected_vec{expected.begin(), expected.end()};
  dealii::Vector<double> result_vec{result.begin(), result.end()};
  return CompareVector(expected_vec, result_vec, tol);
}

AssertionResult CompareMPIMatrices(
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


} // namespace testing

} // namespace bart