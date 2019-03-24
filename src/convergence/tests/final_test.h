#ifndef BART_SRC_CONVERGENCE_MOMENTS_TESTS_H_
#define BART_SRC_CONVERGENCE_MOMENTS_TESTS_H_

#include "convergence/final.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace testing {

// -- CUSTOM ASSERTION FOR COMPARING STATUS STRUCTS --

::testing::AssertionResult CompareStatus(const Status lhs, const Status rhs) {
  if (lhs.iteration_number != rhs.iteration_number) {
    return ::testing::AssertionFailure() << "Iteration number mismatch "
                                         << lhs.iteration_number << " != "
                                         << rhs.iteration_number;
  } else if (lhs.max_iterations != rhs.max_iterations) {
    return ::testing::AssertionFailure() << "Max iteration number mismatch "
                                         << lhs.max_iterations << " != "
                                         << rhs.max_iterations;
  } else if (lhs.is_complete != rhs.is_complete) {
    return ::testing::AssertionFailure() << "Convergence completion mismatch "
                                         << lhs.is_complete << " != "
                                         << rhs.is_complete;
  } else if (lhs.failed_index != rhs.failed_index) {
    return ::testing::AssertionFailure() << "Failed index mismatch "
                                         << lhs.failed_index.value_or(-1) << " != "
                                         << rhs.failed_index.value_or(-1);
  } else if (lhs.delta != rhs.delta) {
    return ::testing::AssertionFailure() << "Delta mismatch "
                                         << lhs.delta.value_or(-1) << " != "
                                         << rhs.delta.value_or(-1);
  }
  return ::testing::AssertionSuccess();
}

template <typename CompareType>
class ConvergenceFinalTest : public ::testing::Test {
  protected:
  void TestBaseMethods(bart::convergence::Final<CompareType>* test_checker) {
    EXPECT_EQ(test_checker->max_iterations(), 100);
    test_checker->SetMaxIterations(50);
    EXPECT_EQ(test_checker->max_iterations(), 50);

    EXPECT_EQ(test_checker->iteration(), 0);
    test_checker->SetIteration(10);
    EXPECT_EQ(test_checker->iteration(), 10);

    // -- ERRORS --
    EXPECT_ANY_THROW(test_checker->SetIteration(-10));
    EXPECT_ANY_THROW(test_checker->SetMaxIterations(-10));
    EXPECT_ANY_THROW(test_checker->SetMaxIterations(0));
  };
};

} // namespace testing

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_TESTS_H_