#ifndef BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_

#include "convergence/final_i.h"
#include "convergence/status.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

template <typename CompareType>
class FinalCheckerMock : public FinalI<CompareType> {
 public:
  using typename FinalI<CompareType>::IterationNumber;
  MOCK_METHOD2_T(CheckFinalConvergence, Status(CompareType& current_iteration,
                                               CompareType& previous_iteration));
  MOCK_CONST_METHOD0_T(convergence_status, Status());
  MOCK_CONST_METHOD0_T(convergence_is_complete, bool());
  MOCK_CONST_METHOD0_T(max_iterations, IterationNumber());
  MOCK_CONST_METHOD0_T(iteration, IterationNumber());
  MOCK_METHOD1_T(SetMaxIterations, FinalCheckerMock&(IterationNumber to_set));
  MOCK_METHOD1_T(SetIteration, FinalCheckerMock&(IterationNumber to_set));
};

} // namespace convergence

} // namespace bart

#endif //BART_SRC_CONVERGENCE_TESTS_FINAL_CHECKER_MOCK_H_
