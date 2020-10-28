#ifndef BART_SRC_SOLVER_GROUP_FACTORY_HPP_
#define BART_SRC_SOLVER_GROUP_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.h"

namespace bart::solver::group {

class SingleGroupSolverI;

enum class GroupSolverName {
  kDefaultImplementation = 0, // solver::group::SingleGroupSolver
};

BART_INTERFACE_FACTORY(SingleGroupSolverI, GroupSolverName)

[[nodiscard]] inline auto to_string(GroupSolverName to_convert) -> std::string {
  switch (to_convert) {
    case GroupSolverName::kDefaultImplementation:
      return std::string{"GroupSolverName::kDefaultImplementation"};
  }
  return std::string{"Unknown GroupSolverName to string conversion requested."};
}


} // namespace bart::solver::group

#endif //BART_SRC_SOLVER_GROUP_FACTORY_HPP_
