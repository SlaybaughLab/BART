#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_

#include <vector>

#include "utility/has_description.h"
#include "iteration/subroutine/subroutine_i.hpp"

namespace bart {

namespace system {
class System;
} // namespace system
namespace iteration {

namespace outer {

class OuterIterationI : public utility::HasDescription {
 public:
  using Subroutine = subroutine::SubroutineI;
  virtual ~OuterIterationI() = default;
  virtual auto IterateToConvergence(system::System &system) -> void = 0;
  virtual auto AddPostIterationSubroutine(std::unique_ptr<Subroutine>) -> OuterIterationI& = 0;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_I_HPP_
