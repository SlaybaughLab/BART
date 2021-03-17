#ifndef BART_SRC_ITERATION_GROUP_TESTS_GROUP_SOLVE_ITERATION_MOCK_HPP_
#define BART_SRC_ITERATION_GROUP_TESTS_GROUP_SOLVE_ITERATION_MOCK_HPP_

#include "iteration/group/group_solve_iteration_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace iteration {

namespace group {

class GroupSolveIterationMock : public GroupSolveIterationI {
 public:
  MOCK_METHOD(void, Iterate, (system::System &system), (override));
  MOCK_METHOD(GroupSolveIterationMock&, UpdateThisAngularSolutionMap,
              (system::solution::EnergyGroupToAngularSolutionPtrMap&), (override));
};

} // namespace group

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_GROUP_TESTS_GROUP_SOLVE_ITERATION_MOCK_HPP_
