#include "iteration/outer/outer_power_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

OuterPowerIteration::OuterPowerIteration(
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    const std::shared_ptr<SourceUpdater> &source_updater_ptr)
    : OuterIteration(
        std::move(convergence_checker_ptr),
        source_updater_ptr) {}
convergence::Status OuterPowerIteration::CheckConvergence(system::System &system) {
  double k_effective_last = 0.0;
  double k_effective_current = 0.0;
  return convergence_checker_ptr_->CheckFinalConvergence(
      k_effective_current, k_effective_last);
}
void OuterPowerIteration::UpdateSystem(system::System &system, const int group,
    const int angle) {
  source_updater_ptr_->UpdateFissionSource(system, group, angle);
}

} // namespace outer

} // namespace iteration

} // namespace bart