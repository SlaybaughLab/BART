#ifndef BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_

#include "convergence/final_i.h"
#include "convergence/status.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

template <typename CompareType>
class FinalCheckerMock : public FinalI<CompareType> {
 public:
  using typename FinalI<CompareType>::IterationNumber;
  MOCK_METHOD(Status, CheckFinalConvergence, (CompareType& current_iteration,
      CompareType& previous_iteration), (override));
  MOCK_METHOD(Status, convergence_status, (), (override, const));
  MOCK_METHOD(bool, convergence_is_complete, (), (override, const));
  MOCK_METHOD(IterationNumber, max_iterations, (), (override, const));
  MOCK_METHOD(IterationNumber, iteration, (), (override, const));
  MOCK_METHOD(FinalCheckerMock&, SetMaxIterations, (IterationNumber to_set),
              (override));
  MOCK_METHOD(FinalCheckerMock&, SetIteration, (IterationNumber to_set),
              (override));
  MOCK_METHOD(void, Reset, (), (override));
};

} // namespace convergence

} // namespace bart

#endif //BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_
