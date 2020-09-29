#ifndef BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_
#define BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_

#include "solver/solver_names.h"
#include "solver/factory/solver_factory.h"

namespace bart {

namespace solver {

namespace factory {

BART_INTERFACE_FACTORY_REGISTRAR(LinearI, LinearSolverName)

} // namespace factory

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_
