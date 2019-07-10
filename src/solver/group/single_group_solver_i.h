#ifndef BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_
#define BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_

namespace bart {

namespace solver {

namespace group {

/*! \brief Interface for solver that solves a single groups.
 *
 * The purpose of these classes is to solve a system for a specific group, using
 * a provided system of equations and given a solution vector to update.
 *
 */
class SingleGroupSolverI {
 public:
  virtual ~SingleGroupSolverI() = default;
};

} // namespace group

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_
