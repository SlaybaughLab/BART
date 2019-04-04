#ifndef BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_
#define BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_

#include <deal.II/lac/vector.h>
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace testing {

using ::testing::AssertionResult;
using ::testing::AssertionFailure;
using ::testing::AssertionSuccess;

AssertionResult CompareVector(const dealii::Vector<double>& expected,
                              const dealii::Vector<double>& result,
                              const double tol = 1e-6) {
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
                              const double tol = 1e-6) {
  dealii::Vector<double> expected_vec{expected.begin(), expected.end()};
  dealii::Vector<double> result_vec{result.begin(), result.end()};
   return CompareVector(expected_vec, result_vec, tol);
}

} // namespace testing

} // namespace bart

#endif // BART_SRC_TEST_HELPERS_TEST_ASSERTIONS_H_