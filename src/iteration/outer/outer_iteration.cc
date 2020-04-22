#include "iteration/outer/outer_iteration.h"

#include "convergence/status.h"

namespace bart {

namespace iteration {

namespace outer {

template <typename ConvergenceType>
OuterIteration<ConvergenceType>::OuterIteration(
    std::unique_ptr<GroupIterator> group_iterator_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    const std::shared_ptr<Reporter> &reporter_ptr)
    : group_iterator_ptr_(std::move(group_iterator_ptr)),
      convergence_checker_ptr_(std::move(convergence_checker_ptr)),
      reporter_ptr_(reporter_ptr) {

  AssertThrow(group_iterator_ptr_ != nullptr,
              dealii::ExcMessage("GroupSolveIteration pointer passed to "
                                 "OuterIteration constructor is null"))

  AssertThrow(convergence_checker_ptr_ != nullptr,
              dealii::ExcMessage("Convergence checker pointer passed to "
                                 "OuterIteration constructor is null"))
}

template <typename ConvergenceType>
void OuterIteration<ConvergenceType>::IterateToConvergence(
    system::System &system) {
  const int total_groups = system.total_groups;
  const int total_angles = system.total_angles;

  convergence::Status convergence_status;

  do {

    if (!convergence_status.is_complete) {
      for (int group = 0; group < total_groups; ++group) {
        for (int angle = 0; angle < total_angles; ++angle) {
          UpdateSystem(system, group, angle);
        }
      }
    }

    InnerIterationToConvergence(system);

    convergence_status = CheckConvergence(system);

    if (reporter_ptr_ != nullptr) {
      reporter_ptr_->Report("Outer iteration Status: ");
      reporter_ptr_->Report(convergence_status);
    }

  } while (!convergence_status.is_complete);
}

template <typename ConvergenceType>
void OuterIteration<ConvergenceType>::InnerIterationToConvergence(
    system::System &system) {
  group_iterator_ptr_->Iterate(system);
}

template class OuterIteration<double>;

} // namespace outer

} // namespace iteration

} // namespace bart