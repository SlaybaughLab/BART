#include "solver/group/single_group_solver.h"

namespace bart {

namespace solver {

namespace group {

SingleGroupSolver::SingleGroupSolver(
    std::unique_ptr<LinearSolver> linear_solver_ptr)
    : linear_solver_ptr_(std::move(linear_solver_ptr)) {}

void SingleGroupSolver::SolveGroup(const int group,
                                   const system::System &system,
                                   system::solution::MPIAngularI &group_solution) {

}

} // namespace group

} // namespace solver

} //namespace bart
