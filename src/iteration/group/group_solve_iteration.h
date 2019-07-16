#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_

#include "convergence/final_i.h"
#include "iteration/group/group_solve_iteration_i.h"
#include "iteration/updater/source_updater_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"

#include <memory>

#include "solver/group/single_group_solver_i.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
class GroupSolveIteration : public GroupSolveIterationI {
 public:
  using GroupSolver = solver::group::SingleGroupSolverI;
  using ConvergenceChecker = convergence::FinalI<system::moments::MomentVector>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsI<dim>;
  using GroupSolution = system::solution::MPIGroupAngularSolutionI;
  using SourceUpdater = iteration::updater::SourceUpdaterI;

  GroupSolveIteration(
      std::unique_ptr<GroupSolver> group_solver_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      std::unique_ptr<MomentCalculator> moment_calculator_ptr,
      std::shared_ptr<GroupSolution> group_solution_ptr,
      std::unique_ptr<SourceUpdater> source_updater_ptr);
  virtual ~GroupSolveIteration() = default;

  void Iterate(system::System &system) override;

  GroupSolver* group_solver_ptr() const {
    return group_solver_ptr_.get();
  }

  ConvergenceChecker* convergence_checker_ptr() const {
    return convergence_checker_ptr_.get();
  }

  MomentCalculator* moment_calculator_ptr() const {
    return moment_calculator_ptr_.get();
  }

  std::shared_ptr<GroupSolution> group_solution_ptr() const {
    return group_solution_ptr_;
  }

  SourceUpdater* source_updater_ptr() const {
    return source_updater_ptr_.get();
  }

 protected:

  virtual void SolveGroup(const int group, system::System &system);
  virtual system::moments::MomentVector GetScalarFlux(const int group,
                                                      system::System& system);
  virtual convergence::Status CheckConvergence(
      system::moments::MomentVector& current_iteration,
      system::moments::MomentVector& previous_iteration);
  virtual void UpdateSystem(system::System& system, const int group,
                            const int angle) = 0;
  virtual void UpdateCurrentMoments(system::System &system, const int group);

  std::unique_ptr<GroupSolver> group_solver_ptr_ = nullptr;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_ = nullptr;
  std::unique_ptr<MomentCalculator> moment_calculator_ptr_ = nullptr;
  std::shared_ptr<GroupSolution> group_solution_ptr_ = nullptr;
  std::unique_ptr<SourceUpdater> source_updater_ptr_ = nullptr;
};

} // namespace group

} // namespace iteration



} //namespace bart

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
