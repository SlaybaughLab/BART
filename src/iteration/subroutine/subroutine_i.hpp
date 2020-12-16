#ifndef BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_

namespace bart::iteration::subroutine {

namespace system {
struct System;
} // namespace system

class SubroutineI {
 public:
  virtual ~SubroutineI() = default;
  virtual auto Execute(system::System&) -> void;
};

} // namespace bart::iteration::subroutine

#endif //BART_SRC_ITERATION_SUBROUTINE_SUBROUTINE_I_HPP_
