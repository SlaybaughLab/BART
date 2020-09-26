#ifndef BART_SRC_SOLVER_GROUP_SINGLE_GROUP_SOLVER_H_
#define BART_SRC_SOLVER_GROUP_SINGLE_GROUP_SOLVER_H_

#include <memory>

#include "solver/group/single_group_solver_i.h"
#include "solver/linear/linear_i.h"

namespace bart {

namespace solver {

namespace group {

class SingleGroupSolver : public SingleGroupSolverI {
 public:

  using LinearSolver = solver::LinearI;

  SingleGroupSolver(std::unique_ptr<LinearSolver> linear_solver_ptr);
  virtual ~SingleGroupSolver() = default;

  void SolveGroup(const int group,
                  const system::System &system,
                  system::solution::MPIGroupAngularSolutionI &group_solution) override;

  LinearSolver* linear_solver_ptr() const {
    return linear_solver_ptr_.get();
  }
 protected:
  std::unique_ptr<LinearSolver> linear_solver_ptr_ = nullptr;

};

} // namespace group

} // namespace solver

} //namespace bart

#endif //BART_SRC_SOLVER_GROUP_SINGLE_GROUP_SOLVER_H_
