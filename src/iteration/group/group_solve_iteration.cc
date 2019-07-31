#include "iteration/group/group_solve_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
GroupSolveIteration<dim>::GroupSolveIteration(
    std::unique_ptr<GroupSolver> group_solver_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<MomentCalculator> moment_calculator_ptr,
    const std::shared_ptr<GroupSolution> &group_solution_ptr,
    const std::shared_ptr<SourceUpdater> &source_updater_ptr,
    const std::shared_ptr<Reporter> &reporter_ptr)
    : group_solver_ptr_(std::move(group_solver_ptr)),
      convergence_checker_ptr_(std::move(convergence_checker_ptr)),
      moment_calculator_ptr_(std::move(moment_calculator_ptr)),
      group_solution_ptr_(group_solution_ptr),
      source_updater_ptr_(source_updater_ptr),
      reporter_ptr_(reporter_ptr) {

  AssertThrow(group_solver_ptr_ != nullptr,
              dealii::ExcMessage("Group solver pointer passed to "
                                 "GroupSolveIteration constructor is null"));
  AssertThrow(convergence_checker_ptr_ != nullptr,
              dealii::ExcMessage("Convergence checker pointer passed to "
                                 "GroupSolveIteration constructor is null"));
  AssertThrow(moment_calculator_ptr_ != nullptr,
              dealii::ExcMessage("Moment calculator pointer passed to "
                                 "GroupSolveIteration constructor is null"));
  AssertThrow(group_solution_ptr_ != nullptr,
              dealii::ExcMessage("Group solution pointer passed to "
                                 "GroupSolveIteration constructor is null"));
  AssertThrow(source_updater_ptr_ != nullptr,
              dealii::ExcMessage("Source updater pointer passed to "
                                 "GroupSolveIteration constructor is null"));


}

template<int dim>
void GroupSolveIteration<dim>::Iterate(system::System &system) {

  const int total_groups = system.current_moments->total_groups();
  const int total_angles = group_solution_ptr_->total_angles();
  system::moments::MomentVector current_scalar_flux, previous_scalar_flux;

  if (reporter_ptr_ != nullptr)
    reporter_ptr_->Report("..Inner group iteration\n");

  for (int group = 0; group < total_groups; ++group) {
    if (reporter_ptr_ != nullptr) {
      std::string report{"....Group: "};
      report += std::to_string(group);
      report += "\n";
      reporter_ptr_->Report(report);
    }

    convergence::Status convergence_status;
    convergence_checker_ptr_->Reset();
    do {

      if (!convergence_status.is_complete) {
        for (int angle = 0; angle < total_angles; ++angle)
          UpdateSystem(system, group, angle);
      }

      previous_scalar_flux = current_scalar_flux;

      SolveGroup(group, system);

      current_scalar_flux = GetScalarFlux(group, system);

      if (convergence_status.iteration_number == 0) {
        previous_scalar_flux = current_scalar_flux;
        previous_scalar_flux = 0;
      }

      convergence_status = CheckConvergence(current_scalar_flux,
                                            previous_scalar_flux);

      if (reporter_ptr_ != nullptr)
        reporter_ptr_->Report(convergence_status);

    } while (!convergence_status.is_complete);
    UpdateCurrentMoments(system, group);
  }
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

template <int dim>
convergence::Status GroupSolveIteration<dim>::CheckConvergence(
    system::moments::MomentVector &current_iteration,
    system::moments::MomentVector &previous_iteration) {
  return convergence_checker_ptr_->CheckFinalConvergence(current_iteration,
                                                         previous_iteration);
}

template <int dim>
void GroupSolveIteration<dim>::UpdateCurrentMoments(system::System &system,
                                                    const int group) {
  auto& current_moments = *system.current_moments;
  const int max_harmonic_l = current_moments.max_harmonic_l();


  for (int l = 0; l <= max_harmonic_l; ++l) {
    for (int m = -l; m <= l; ++m) {
      current_moments[{group, l, m}] = moment_calculator_ptr_->CalculateMoment(
          group_solution_ptr_.get(), group, l, m);
    }
  }
}

template class GroupSolveIteration<1>;
template class GroupSolveIteration<2>;
template class GroupSolveIteration<3>;

} // namespace group

} // namespace iteration

} // namespace bart
