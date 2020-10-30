#ifndef BART_SRC_SOLVER_BUILDER_SOLVER_BUILDER_HPP_
#define BART_SRC_SOLVER_BUILDER_SOLVER_BUILDER_HPP_

#include "solver/group/single_group_solver_i.h"

namespace bart::solver::builder {

enum class SolverName {
  kDefaultGMRESGroupSolver = 0,
};

class SolverBuilder {
 public:
  template <typename ...Args>
  [[nodiscard]] static auto BuildSolver(const SolverName, Args...) -> std::unique_ptr<group::SingleGroupSolverI>;
};

} // namespace bart::solver::builder

#endif //BART_SRC_SOLVER_BUILDER_SOLVER_BUILDER_HPP_
