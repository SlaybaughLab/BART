#ifndef BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_
#define BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_

#include "iteration/subroutine/subroutine_i.hpp"

namespace bart::iteration::subroutine {

class GetScalarFluxFromFramework : public SubroutineI {
 public:
  auto Execute(system::System &) -> void override {};
};

} // namespace bart::iteration::subroutine

#endif //BART_SRC_ITERATION_SUBROUTINE_GET_SCALAR_FLUX_FROM_FRAMEWORK_HPP_
