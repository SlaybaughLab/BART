#include "solver/group/single_group_solver.h"

#include "system/system.h"
#include "system/solution/mpi_angular_i.h"

namespace bart {

namespace solver {

namespace group {

SingleGroupSolver::SingleGroupSolver(
    std::unique_ptr<LinearSolver> linear_solver_ptr)
    : linear_solver_ptr_(std::move(linear_solver_ptr)) {}

void SingleGroupSolver::SolveGroup(const int group,
                                   const system::System &system,
                                   system::solution::MPIAngularI &group_solution) {
  const int total_angles = group_solution.total_angles();
  AssertThrow(total_angles > 0,
      dealii::ExcMessage("Error in SolveGroup, total angles provided by group "
                         "solution must be > 0"));
  AssertThrow(group >= 0,
      dealii::ExcMessage("Error in SolveGroup, invalid group index provided, "
                         "value is less than zero"));
  AssertThrow(group < group_solution.total_groups(),
      dealii::ExcMessage("Error in SolveGroup, invalid group index provided, "
                         "not within range of groups in solution"));
}

} // namespace group

} // namespace solver

} //namespace bart
