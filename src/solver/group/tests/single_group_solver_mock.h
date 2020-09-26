#ifndef BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_MOCK_H_
#define BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_MOCK_H_

#include "solver/group/single_group_solver_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace solver {

namespace group {

class SingleGroupSolverMock : public SingleGroupSolverI {
 public:
  MOCK_METHOD(void, SolveGroup,
              (const int group,
                  const system::System& system,
                  system::solution::MPIGroupAngularSolutionI& group_solution),
              (override));
};

} // namespace group

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_MOCK_H_
