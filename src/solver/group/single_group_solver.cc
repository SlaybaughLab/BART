#include "solver/group/single_group_solver.h"

#include "solver/group/factory.hpp"
#include "system/system.h"
#include "system/solution/mpi_group_angular_solution_i.h"

namespace bart {

namespace solver {

namespace group {

SingleGroupSolver::SingleGroupSolver(
    std::unique_ptr<LinearSolver> linear_solver_ptr)
    : linear_solver_ptr_(std::move(linear_solver_ptr)) {}

bool SingleGroupSolver::is_registered_ =
    SingleGroupSolverIFactory<std::unique_ptr<LinearSolver>>::get()
    .RegisterConstructor(GroupSolverName::kDefaultImplementation,
        [](std::unique_ptr<LinearSolver> linear_solver_ptr) {
          std::unique_ptr<SingleGroupSolverI> return_ptr =
              std::make_unique<SingleGroupSolver>(std::move(linear_solver_ptr));
          return return_ptr; });

void SingleGroupSolver::SolveGroup(const int group,
                                   const system::System &system,
                                   system::solution::MPIGroupAngularSolutionI &group_solution) {
  const int total_angles = group_solution.total_angles();
  AssertThrow(total_angles > 0,
      dealii::ExcMessage("Error in SolveGroup, total angles provided by group "
                         "solution must be > 0"));
  AssertThrow(group >= 0,
      dealii::ExcMessage("Error in SolveGroup, invalid group index provided, "
                         "value is less than zero"));

  for (int angle = 0; angle < total_angles; ++angle) {
    system::Index index{group, angle};
    auto& solution = group_solution[angle];
    auto left_hand_side_ptr = system.left_hand_side_ptr_->GetFullTermPtr(index);
    auto right_hand_side_ptr = system.right_hand_side_ptr_->GetFullTermPtr(index);
    dealii::PETScWrappers::PreconditionNone no_conditioner(*left_hand_side_ptr);

    linear_solver_ptr_->Solve(
        left_hand_side_ptr.get(),
        &solution,
        right_hand_side_ptr.get(),
        &no_conditioner);
  }
}

} // namespace group

} // namespace solver

} //namespace bart
