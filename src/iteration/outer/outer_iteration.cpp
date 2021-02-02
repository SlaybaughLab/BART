#include "iteration/outer/outer_iteration.hpp"

#include "convergence/status.hpp"

namespace bart::iteration::outer {

template <typename ConvergenceType>
OuterIteration<ConvergenceType>::OuterIteration(std::unique_ptr<GroupIterator> group_iterator_ptr,
                                                std::unique_ptr<ConvergenceChecker> convergence_checker_ptr)
    : group_iterator_ptr_(std::move(group_iterator_ptr)),
      convergence_checker_ptr_(std::move(convergence_checker_ptr)) {

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
  bool is_complete{ false };
  do {
    is_complete = Iterate(system);
    if (post_iteration_subroutine_ptr_ != nullptr)
      post_iteration_subroutine_ptr_->Execute(system);
  } while (!is_complete);
}

template <typename ConvergenceType>
void OuterIteration<ConvergenceType>::InnerIterationToConvergence(
    system::System &system) {
  group_iterator_ptr_->Iterate(system);
}

template<typename ConvergenceType>
auto OuterIteration<ConvergenceType>::Iterate(system::System &system) -> bool {
  for (int group = 0; group < system.total_groups; ++group) {
    for (int angle = 0; angle < system.total_angles; ++angle) {
      UpdateSystem(system, group, angle);
    }
  }

  InnerIterationToConvergence(system);

  auto convergence_status = CheckConvergence(system);
  if (convergence_status.delta.has_value()) {
    data_names::IterationErrorPort::Expose({convergence_status.iteration_number,
                                            convergence_status.delta.value()});
  }

  data_names::StatusPort::Expose("Outer iteration Status: ");
  data_names::ConvergenceStatusPort::Expose(convergence_status);
  data_names::SolutionMomentsPort::Expose(*system.current_moments);

  return convergence_status.is_complete;
}

template class OuterIteration<double>;

} // namespace bart::iteration::outer