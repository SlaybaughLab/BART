#include "iteration/group/group_source_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
GroupSourceIteration<dim>::GroupSourceIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    std::shared_ptr<GroupSolution> group_solution_ptr,
    std::unique_ptr<SourceUpdater> source_updater_ptr)
    : GroupSolveIteration<dim>(std::move(group_solver_ptr),
        std::move(convergence_checker_ptr),
        std::move(moment_calculator_ptr),
        group_solution_ptr,
        std::move(source_updater_ptr))
{}
template<int dim>
void GroupSourceIteration<dim>::UpdateSystem(system::System &system, const int group) {

}

template class GroupSourceIteration<1>;
template class GroupSourceIteration<2>;
template class GroupSourceIteration<3>;


} // namespace group

} // namespace iteration

} // namespace bart
