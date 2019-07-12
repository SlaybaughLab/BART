#include "iteration/group/group_solve_iteration.h"

namespace bart {

namespace iteration {

namespace group {

GroupSolveIteration::GroupSolveIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr)
    : group_solver_ptr_(std::move(group_solver_ptr))
    {}

} // namespace group

} // namespace iteration

} // namespace bart
