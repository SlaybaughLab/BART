#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_

#include <memory>

#include "iteration/outer/outer_iteration_i.h"
#include "iteration/updater/source_updater_i.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterIteration : public OuterIterationI {
 public:
  using SourceUpdater = iteration::updater::SourceUpdaterI;

  OuterIteration(
      const std::shared_ptr<SourceUpdater> &source_updater_ptr);
  virtual ~OuterIteration() = default;
  virtual void IterateToConvergence(system::System &system) {}

  SourceUpdater* source_updater_ptr() const {
    return source_updater_ptr_.get();
  }

 protected:
  std::shared_ptr<SourceUpdater> source_updater_ptr_ = nullptr;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
