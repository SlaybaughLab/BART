#ifndef BART_SRC_SOLVER_GROUP_FACTORY_HPP_
#define BART_SRC_SOLVER_GROUP_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.h"

namespace bart::solver::group {

class SingleGroupSolverI;

enum class GroupSolverName {
  kDefaultImplementation = 0, // solver::group::SingleGroupSolver
};

BART_INTERFACE_FACTORY(SingleGroupSolverI, GroupSolverName)

} // namespace bart::solver::group

#endif //BART_SRC_SOLVER_GROUP_FACTORY_HPP_
