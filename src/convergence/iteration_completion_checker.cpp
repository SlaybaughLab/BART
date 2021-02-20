#include "convergence/iteration_completion_checker.hpp"
#include "system/moments/spherical_harmonic_types.h"

namespace bart::convergence {

template <typename CompareType>
IterationCompletionChecker<CompareType>::IterationCompletionChecker(
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr)
    : convergence_checker_ptr_(std::move(convergence_checker_ptr)) {}

template<typename CompareType>
auto IterationCompletionChecker<CompareType>::ConvergenceStatus(const CompareType &current_iteration,
                                                                const CompareType &previous_iteration) -> Status {
  // Get convergence status
  convergence_status_.is_complete = convergence_checker_ptr_->IsConverged(current_iteration, previous_iteration);
  convergence_status_.delta = convergence_checker_ptr_->delta();
  if (!convergence_status_.is_complete) {
    convergence_status_.failed_index = convergence_checker_ptr_->failed_index();
  } else {
    convergence_status_.failed_index = std::nullopt;
  }
  // Increment iteration number
  ++convergence_status_.iteration_number;
  // Check if max iterations have been reached
  if (convergence_status_.iteration_number == convergence_status_.max_iterations)
    convergence_status_.is_complete = true;
  return convergence_status_;
}

template<typename CompareType>
auto IterationCompletionChecker<CompareType>::SetMaxIterations(IterationNumber to_set)
-> IterationCompletionChecker<CompareType>& {
  AssertThrow(to_set > 0, dealii::ExcMessage("Max iterations must be > 0"))
  convergence_status_.max_iterations = to_set;
  return *this;
}

template<typename CompareType>
auto IterationCompletionChecker<CompareType>::SetIteration(IterationNumber to_set)
-> IterationCompletionChecker<CompareType> & {
  AssertThrow(to_set >= 0, dealii::ExcMessage("Iteration must be >= 0"));
  convergence_status_.iteration_number = to_set;
  return *this;
}

template<typename CompareType>
auto IterationCompletionChecker<CompareType>::Reset() -> void {
  Status convergence_status;
  convergence_status.max_iterations = convergence_status_.max_iterations;
  convergence_status_ = convergence_status;
}

template class IterationCompletionChecker<system::moments::MomentVector>;
template class IterationCompletionChecker<system::moments::MomentsMap>;
template class IterationCompletionChecker<double>;

} // namespace bart::convergence
