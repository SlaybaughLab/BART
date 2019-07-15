#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOURCE_ITERATION_H_
#define BART_SRC_ITERATION_GROUP_GROUP_SOURCE_ITERATION_H_

#include "iteration/group/group_solve_iteration.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
class GroupSourceIteration : public GroupSolveIteration<dim> {
 public:
  using typename GroupSolveIteration<dim>::GroupSolver;
  using typename GroupSolveIteration<dim>::ConvergenceChecker;
  using typename GroupSolveIteration<dim>::MomentCalculator;
  using typename GroupSolveIteration<dim>::GroupSolution;
  using typename GroupSolveIteration<dim>::SourceUpdater;

  GroupSourceIteration(
      std::unique_ptr<GroupSolver> group_solver_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      std::unique_ptr<MomentCalculator> moment_calculator_ptr,
      std::shared_ptr<GroupSolution> group_solution_ptr,
      std::unique_ptr<SourceUpdater> source_updater_ptr);
  virtual ~GroupSourceIteration() = default;

 protected:
  void UpdateSystem(system::System &system, const int group) override;

};

} // namespace group

} // namespace iteration

} //namespace bart

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOURCE_ITERATION_H_
