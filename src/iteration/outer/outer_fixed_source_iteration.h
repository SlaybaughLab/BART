#ifndef BART_SRC_ITERATION_OUTER_OUTER_FIXED_SOURCE_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_FIXED_SOURCE_ITERATION_H_

#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterFixedSourceIteration : public OuterIteration<double> {
 public:
  using typename OuterIteration<double>::GroupIterator;
  using typename OuterIteration<double>::ConvergenceChecker;

  OuterFixedSourceIteration(
      std::unique_ptr<GroupIterator> group_iterator_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr);
  virtual ~OuterFixedSourceIteration() = default;
  convergence::Status CheckConvergence(system::System &system) override;
  void UpdateSystem(system::System &/*system*/,
                    const int /*group*/,
                    const int /*angle*/) override {}
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_FIXED_SOURCE_ITERATION_H_
