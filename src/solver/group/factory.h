#ifndef BART_SRC_SOLVER_GROUP_FACTORY_H_
#define BART_SRC_SOLVER_GROUP_FACTORY_H_

#include "utility/factory/auto_registering_factory.h"

namespace bart {

namespace solver {

namespace group {

class SingleGroupSolverI;

enum class GroupSolverName {
  kDefaultImplementation = 0, // solver::group::SingleGroupSolver
};

BART_INTERFACE_FACTORY(SingleGroupSolverI, GroupSolverName)

} // namespace group

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_GROUP_FACTORY_H_
