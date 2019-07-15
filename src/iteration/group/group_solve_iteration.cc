#include "iteration/group/group_solve_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
GroupSolveIteration<dim>::GroupSolveIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    std::shared_ptr<GroupSolution> group_solution_ptr,
    std::unique_ptr<SourceUpdater> source_updater_ptr)
    : group_solver_ptr_(std::move(group_solver_ptr)),
      convergence_checker_ptr_(std::move(convergence_checker_ptr)),
      moment_calculator_ptr_(std::move(moment_calculator_ptr)),
      group_solution_ptr_(group_solution_ptr),
      source_updater_ptr_(std::move(source_updater_ptr))
{}
template<int dim>
void GroupSolveIteration<dim>::Iterate(system::System &system) {
  convergence::Status convergence_status;

  const int total_groups = system.current_moments->total_groups();
  system::moments::MomentVector current_scalar_flux, previous_scalar_flux;

  do {
    for (int group = 0; group < total_groups; ++group) {
      SolveGroup(group, system);

      current_scalar_flux = GetScalarFlux(group, system);

      convergence_status.is_complete = true;
    }
  } while (!convergence_status.is_complete);
}

template <int dim>
void GroupSolveIteration<dim>::SolveGroup(int group, system::System &system) {
  group_solver_ptr_->SolveGroup(group, system, *group_solution_ptr_);
}

template <int dim>
system::moments::MomentVector GroupSolveIteration<dim>::GetScalarFlux(
    const int group, system::System &) {
  return moment_calculator_ptr_->CalculateMoment(group_solution_ptr_.get(),
                                                 group, 0, 0);
}

template class GroupSolveIteration<1>;
template class GroupSolveIteration<2>;
template class GroupSolveIteration<3>;

} // namespace group

} // namespace iteration

} // namespace bart
