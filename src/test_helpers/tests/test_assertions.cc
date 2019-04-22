#include "test_helpers/test_assertions.h"

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

  for (unsigned int i = first_local_row; i < last_local_row; ++i) {
    for (unsigned int j = 0; j < n_columns; ++j) {
      if (result(i, j) != expected(i, j)) {
        return AssertionFailure() << "Entry (" << i << ", " << j <<
                                  ") has value: " << result.el(i, j) <<
                                  ", expected: " << expected.el(i, j);
      }
    }
  }
  return AssertionSuccess();
}


} // namespace testing

} // namespace bart