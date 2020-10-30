#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_

#include <vector>

#include "utility/has_description.h"

namespace bart {

namespace system {
class System;
} // namespace system
namespace iteration {

namespace outer {

class OuterIterationI : public utility::HasDescription {
 public:
  virtual ~OuterIterationI() = default;
  virtual void IterateToConvergence(system::System &system) = 0;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_
