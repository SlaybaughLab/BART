#ifndef BART_SRC_SOLVER_SOLVER_NAMES_H_
#define BART_SRC_SOLVER_SOLVER_NAMES_H_

namespace bart {

namespace solver {

// Class names for classes in the group namespace that solve a single group
enum class GroupSolverName {
  kDefaultImplementation = 0, // solver::group::SingleGroupSolver
};

enum class LinearSolverName {
  kGMRES = 0, //solver::linear::GMRES
};

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_SOLVER_NAMES_H_
