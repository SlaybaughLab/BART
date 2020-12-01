#ifndef BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
#define BART_SRC_SOLVER_LINEAR_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"
#include "solver/linear/linear_i.hpp"

namespace bart::solver::linear {

class LinearI;

enum class LinearSolverName {
  kGMRES = 0, //solver::linear::GMRES
};

BART_INTERFACE_FACTORY(LinearI, LinearSolverName)

// LCOV_EXCL_START
[[nodiscard]] inline auto to_string(LinearSolverName to_convert) -> std::string {
  switch (to_convert) {
    case LinearSolverName::kGMRES:
      return std::string{"LinearSolverName::kGMRES"};
  }
  return std::string{"String not defined for specified LinearSolverName"};
}
// LCOV_EXCL_STOP

} // namespace bart::solver::linear

#endif //BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
