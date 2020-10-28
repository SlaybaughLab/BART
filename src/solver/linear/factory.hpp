#ifndef BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
#define BART_SRC_SOLVER_LINEAR_FACTORY_HPP_

#include "utility/factory/auto_registering_factory.h"

namespace bart {

namespace solver {

class LinearI;

enum class LinearSolverName {
  kGMRES = 0, //solver::linear::GMRES
};

BART_INTERFACE_FACTORY(LinearI, LinearSolverName)

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_LINEAR_FACTORY_HPP_
