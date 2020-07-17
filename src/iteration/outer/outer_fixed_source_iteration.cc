#include "iteration/outer/outer_fixed_source_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

OuterFixedSourceIteration::OuterFixedSourceIteration(
    std::unique_ptr<GroupIterator> group_iterator_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    const std::shared_ptr<Reporter>& reporter_ptr)
    : OuterIteration<double>(std::move(group_iterator_ptr),
                             std::move(convergence_checker_ptr),
                             reporter_ptr) {}

convergence::Status OuterFixedSourceIteration::CheckConvergence(system::System &system) {
  convergence::Status return_status;
  return_status.is_complete = true;
  return_status.delta = std::nullopt;
  return return_status;
}

} // namespace out

} // namespace iteration

} // namespace bart
