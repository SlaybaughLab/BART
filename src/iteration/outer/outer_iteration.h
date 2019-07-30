#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_

#include "iteration/outer/outer_iteration_i.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterIteration : public OuterIterationI {
 public:
  virtual ~OuterIteration() = default;
  virtual void IterateToConvergence(system::System &system) {}
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
