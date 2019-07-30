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

} // namespace outer

} // namespace iteration

} // namespace bart