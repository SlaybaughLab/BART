#ifndef BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_
#define BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_

#include "system/system.h"
#include "system/solution/mpi_group_angular_solution_i.h"

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
  virtual void SolveGroup(const int group,
                          const system::System& system,
                          system::solution::MPIGroupAngularSolutionI& group_solution) = 0;
};

} // namespace group

} // namespace solver

} // namespace bart

#endif //BART_SRC_SOLVER_GROUP_TESTS_SINGLE_GROUP_SOLVER_I_H_
