#include "iteration/outer/outer_iteration.h"

#include "convergence/status.h"

namespace bart {

namespace iteration {

namespace outer {

template <typename ConvergenceType>
OuterIteration<ConvergenceType>::OuterIteration(
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    const std::shared_ptr<SourceUpdater> &source_updater_ptr)
    : convergence_checker_ptr_(std::move(convergence_checker_ptr)),
      source_updater_ptr_(source_updater_ptr) {

  AssertThrow(convergence_checker_ptr_ != nullptr,
              dealii::ExcMessage("Convergence checker pointer passed to "
                                 "OuterIteration constructor is null"));

  AssertThrow(source_updater_ptr_ != nullptr,
      dealii::ExcMessage("Source updater pointer passed to OuterIteration "
                         "constructor is null"));
}

template <typename ConvergenceType>
void OuterIteration<ConvergenceType>::IterateToConvergence(
    system::System &system) {
  const int total_groups = system.total_groups;
  const int total_angles = system.total_angles;

  convergence::Status convergence_status;

  do {
    convergence_status = CheckConvergence(system);
  } while (!convergence_status.is_complete);
}

template class OuterIteration<double>;

} // namespace outer

} // namespace iteration

} // namespace bart