#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

OuterIteration::OuterIteration(
    const std::shared_ptr<SourceUpdater> &source_updater_ptr)
    : source_updater_ptr_(source_updater_ptr) {
  AssertThrow(source_updater_ptr_ != nullptr,
      dealii::ExcMessage("Source updater pointer passed to OuterIteration "
                         "constructor is null"));
}

} // namespace outer

} // namespace iteration

} // namespace bart