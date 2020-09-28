#ifndef BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_
#define BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_

#include "solver/solver_names.h"
#include "solver/factory/solver_factory.h"

namespace bart {

namespace solver {

namespace factory {

template <typename T, typename ...U>
class SolverFactoryRegistrar {
 public:
  SolverFactoryRegistrar(const LinearSolverName name,
                         const LinearSolverConstructor<U...>& constructor) {
    SolverFactory<U...>::get().RegisterConstructor(name, constructor);
  }
};

} // namespace factory

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_REGISTRAR_H_
