#ifndef BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_

#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterPowerIteration : public OuterIteration<double> {
 public:
  using ConvergenceChecker = convergence::FinalI<double>;
  using OuterIteration<double>::SourceUpdater;

  OuterPowerIteration(
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      const std::shared_ptr<SourceUpdater> &source_updater_ptr);
  virtual ~OuterPowerIteration() = default;
 protected:
  convergence::Status CheckConvergence(system::System &system) override;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
