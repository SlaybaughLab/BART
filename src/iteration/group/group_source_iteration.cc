#include "iteration/group/group_source_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
GroupSourceIteration<dim>::GroupSourceIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    const std::shared_ptr<GroupSolution> &group_solution_ptr,
    const std::shared_ptr<SourceUpdater> &source_updater_ptr,
    const std::shared_ptr<Reporter> &reporter_ptr)
    : GroupSolveIteration<dim>(std::move(group_solver_ptr),
        std::move(convergence_checker_ptr),
        std::move(moment_calculator_ptr),
        group_solution_ptr,
        reporter_ptr) {
  source_updater_ptr_ = source_updater_ptr;
  AssertThrow(source_updater_ptr_ != nullptr,
              dealii::ExcMessage("Source updater pointer passed to "
                                 "GroupSolveIteration constructor is null"));
}

template<int dim>
void GroupSourceIteration<dim>::UpdateSystem(system::System &system,
    const int group, const int angle) {
  this->source_updater_ptr_->UpdateScatteringSource(system,
      system::EnergyGroup(group),
      quadrature::QuadraturePointIndex(angle));
}

template class GroupSourceIteration<1>;
template class GroupSourceIteration<2>;
template class GroupSourceIteration<3>;


} // namespace group

} // namespace iteration

} // namespace bart
