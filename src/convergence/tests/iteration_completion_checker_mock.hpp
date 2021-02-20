#ifndef BART_SRC_CONVERGENCE_TESTS_ITERATION_COMPLETION_CHECKER_MOCK_HPP_
#define BART_SRC_CONVERGENCE_TESTS_ITERATION_COMPLETION_CHECKER_MOCK_HPP_

#include "convergence/iteration_completion_checker_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::convergence {

template <typename CompareType>
class IterationCompletionCheckerMock : public IterationCompletionCheckerI<CompareType> {
 public:
  using typename IterationCompletionCheckerI<CompareType>::IterationNumber;
  MOCK_METHOD(Status, ConvergenceStatus, (const CompareType&, const CompareType&), (override));
  MOCK_METHOD(Status, convergence_status, (), (const, override));
  MOCK_METHOD(bool, convergence_is_complete, (), (const, override));
  MOCK_METHOD(IterationNumber, max_iterations, (), (const, override));
  MOCK_METHOD(IterationNumber, iteration, (), (const, override));
  MOCK_METHOD(IterationCompletionCheckerMock&, SetMaxIterations, (IterationNumber to_set), (override));
  MOCK_METHOD(IterationCompletionCheckerMock&, SetIteration, (IterationNumber to_set), (override));
  MOCK_METHOD(void, Reset, (), (override));
};

} // namespace bart::convergence

#endif //BART_SRC_CONVERGENCE_TESTS_ITERATION_COMPLETION_CHECKER_MOCK_HPP_
