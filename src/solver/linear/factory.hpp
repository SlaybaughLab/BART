#ifndef BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
#define BART_SRC_SOLVER_LINEAR_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.hpp"

namespace bart::solver {

class LinearI;

enum class LinearSolverName {
  kGMRES = 0, //solver::linear::GMRES
};

BART_INTERFACE_FACTORY(LinearI, LinearSolverName)

[[nodiscard]] inline auto to_string(LinearSolverName to_convert) -> std::string {
  switch (to_convert) {
    case LinearSolverName::kGMRES:
      return std::string{"LinearSolverName::kGMRES"};
  }
  return std::string{"Unknown LinearSolverName requested for conversion."};
}

} // namespace bart::solver

#endif //BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
