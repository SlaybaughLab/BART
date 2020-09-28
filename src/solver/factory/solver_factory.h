#ifndef BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_
#define BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_

#include <memory>
#include <unordered_map>

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
class SolverFactory {
 public:

  static SolverFactory& get() {
    static SolverFactory instance;
    return instance;
  }
  bool RegisterConstructor(
      const LinearSolverName name,
      const LinearSolverConstructor<T...>& constructor) {
    return linear_solver_constructors_.insert(
        std::make_pair(name, constructor)).second;
  }
  LinearSolverConstructor<T...> GetConstructor(
      const LinearSolverName name) {
    return linear_solver_constructors_.at(name); }

 private:
  SolverFactory() = default;
  SolverFactory(const SolverFactory&);
  ~SolverFactory() = default;
  std::unordered_map<LinearSolverName, LinearSolverConstructor<T...>>
      linear_solver_constructors_;
};

} // namespace factory

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_FACTORY_SOLVER_FACTORY_H_
