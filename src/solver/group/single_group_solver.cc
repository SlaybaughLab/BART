#include "solver/group/single_group_solver.h"

namespace bart {

namespace solver {

namespace group {

SingleGroupSolver::SingleGroupSolver(
    std::unique_ptr<LinearSolver> linear_solver_ptr)
    : linear_solver_ptr_(std::move(linear_solver_ptr)) {}

} // namespace group

} // namespace solver

} //namespace bart
