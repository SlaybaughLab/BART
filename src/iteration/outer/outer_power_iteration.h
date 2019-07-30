#ifndef BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_

#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterPowerIteration : public OuterIteration {
 public:
  OuterPowerIteration(
      const std::shared_ptr<SourceUpdater> &source_updater_ptr);
  virtual ~OuterPowerIteration() = default;

};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
