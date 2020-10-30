#include "solver/builder/solver_builder.hpp"

#include "solver/group/factory.hpp"
#include "solver/linear/factory.hpp"

namespace bart::solver::builder {

template <>
auto SolverBuilder::BuildSolver(const SolverName) -> std::unique_ptr<group::SingleGroupSolverI> {
  return nullptr;
}

} // namespace bart::solver::builder
