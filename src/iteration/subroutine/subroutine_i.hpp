#ifndef BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_

#include "system/system.hpp"

namespace bart::iteration::subroutine {

enum class SubroutineName {
  kGetScalarFluxFromFramework = 0,
};

class SubroutineI {
 public:
  virtual ~SubroutineI() = default;
  virtual auto Execute(system::System&) -> void = 0;
};

} // namespace bart::iteration::subroutine

#endif //BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_
