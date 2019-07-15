#include "iteration/group/group_solve_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
GroupSolveIteration<dim>::GroupSolveIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    std::shared_ptr<GroupSolution> group_solution_ptr)
    : group_solver_ptr_(std::move(group_solver_ptr)),
      convergence_checker_ptr_(std::move(convergence_checker_ptr)),
      moment_calculator_ptr_(std::move(moment_calculator_ptr)),
      group_solution_ptr_(group_solution_ptr)
{}

template class GroupSolveIteration<1>;
template class GroupSolveIteration<2>;
template class GroupSolveIteration<3>;

} // namespace group

} // namespace iteration

} // namespace bart
