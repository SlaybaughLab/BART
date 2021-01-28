#include "iteration/group/group_source_iteration.hpp"

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
    std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr)
    : GroupSolveIteration<dim>(std::move(group_solver_ptr),
        std::move(convergence_checker_ptr),
        std::move(moment_calculator_ptr),
        group_solution_ptr,
        std::move(moment_map_convergence_checker_ptr)) {
  source_updater_ptr_ = source_updater_ptr;
  AssertThrow(source_updater_ptr_ != nullptr,
              dealii::ExcMessage("Source updater pointer passed to "
                                 "GroupSolveIteration constructor is null"));
  this->set_description("Group source iteration",
                        utility::DefaultImplementation(true));
}

template<int dim>
GroupSourceIteration<dim>::GroupSourceIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    const std::shared_ptr<GroupSolution> &group_solution_ptr,
    const std::shared_ptr<SourceUpdater> &source_updater_ptr,
    const std::shared_ptr<BoundaryConditionsUpdater> &boundary_condition_updater_ptr,
    std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr)
    : GroupSourceIteration<dim>::GroupSourceIteration(
        std::move(group_solver_ptr),
        std::move(convergence_checker_ptr),
        std::move(moment_calculator_ptr),
        group_solution_ptr,
        source_updater_ptr,
        std::move(moment_map_convergence_checker_ptr)) {
  boundary_condition_updater_ptr_ = boundary_condition_updater_ptr;
  AssertThrow(boundary_condition_updater_ptr_ != nullptr,
      dealii::ExcMessage("Boundary conditions updater pointer passed to "
                         "GroupSolveIteration constructor is null"));
  this->set_description("Group source iteration w/boundary conditions update",
                        utility::DefaultImplementation(true));
}

template<int dim>
void GroupSourceIteration<dim>::UpdateSystem(system::System &system,
    const int group, const int angle) {
  this->source_updater_ptr_->UpdateScatteringSource(system,
      system::EnergyGroup(group),
      quadrature::QuadraturePointIndex(angle));
}
template<int dim>
void GroupSourceIteration<dim>::PerformPerGroup(system::System &system,
                                                const int group) {
  GroupSolveIteration<dim>::PerformPerGroup(system, group);
  if (boundary_condition_updater_ptr_ != nullptr) {
    for (int angle = 0; angle < system.total_angles; ++angle) {
      boundary_condition_updater_ptr_->UpdateBoundaryConditions(
          system,
          system::EnergyGroup(group),
          quadrature::QuadraturePointIndex(angle));
    }
  }
}
template<int dim>
auto GroupSourceIteration<dim>::ExposeIterationData(system::System &system) -> void {
  GroupSolveIteration<dim>::ExposeIterationData(system);
  source_updater_ptr_->Expose(source_updater_ptr_->value());
}

template class GroupSourceIteration<1>;
template class GroupSourceIteration<2>;
template class GroupSourceIteration<3>;


} // namespace group

} // namespace iteration

} // namespace bart
