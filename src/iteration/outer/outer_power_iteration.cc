#include "iteration/outer/outer_power_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

OuterPowerIteration::OuterPowerIteration(
    const std::shared_ptr<SourceUpdater> &source_updater_ptr)
    : OuterIteration(source_updater_ptr) {}

} // namespace outer

} // namespace iteration

} // namespace bart