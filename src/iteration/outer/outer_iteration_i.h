#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_H_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_H_

#include <vector>

namespace bart {

namespace system {
class System;
} // namespace system
namespace iteration {

namespace outer {

class OuterIterationI {
 public:
  virtual ~OuterIterationI() = default;
  virtual void IterateToConvergence(system::System &system) = 0;
  virtual std::vector<double> iteration_error() const = 0;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_H_
