#ifndef BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_
#define BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_

#include <memory>
#include <unordered_map>

#include "utility/factory/auto_registering_factory.h"
#include "solver/solver_names.h"

namespace bart {

namespace solver {

// Forward declaration of interfaces built by this factory
namespace group { class SingleGroupSolverI; } // namespace group
class LinearI;

namespace factory {

template <typename ...T>
using LinearSolverConstructor = std::unique_ptr<LinearI>(*)(T...);

template <typename ...T>
class SolverFactory
    : public utility::factory::AutoRegisteringFactory<
        LinearSolverName,
        LinearSolverConstructor<T...>> {};

} // namespace factory

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_
