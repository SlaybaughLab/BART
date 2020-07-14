#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_

#include "convergence/final_i.h"
#include "convergence/reporter/mpi_i.h"
#include "iteration/group/group_solve_iteration_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"

#include <memory>

#include "solver/group/single_group_solver_i.h"
#include "system/solution/solution_types.h"

namespace bart {

namespace iteration {

namespace group {

template <int dim>
class GroupSolveIteration : public GroupSolveIterationI {
 public:
  using GroupSolver = solver::group::SingleGroupSolverI;
  using ConvergenceChecker = convergence::FinalI<system::moments::MomentVector>;
  using MomentMapConvergenceChecker = convergence::FinalI<system::moments::MomentsMap>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsI;
  using GroupSolution = system::solution::MPIGroupAngularSolutionI;
  using Reporter = convergence::reporter::MpiI;
  using EnergyGroupToAngularSolutionPtrMap = system::solution::EnergyGroupToAngularSolutionPtrMap;

  GroupSolveIteration(
      std::unique_ptr<GroupSolver> group_solver_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      std::unique_ptr<MomentCalculator> moment_calculator_ptr,
      const std::shared_ptr<GroupSolution> &group_solution_ptr,
      const std::shared_ptr<Reporter> &reporter_ptr = nullptr,
      std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr = nullptr);

  GroupSolveIteration& UpdateThisAngularSolutionMap(
      EnergyGroupToAngularSolutionPtrMap& to_update) {
    is_storing_angular_solution_ = true;
    angular_solution_ptr_map_ = to_update;
    return *this;
  }

  virtual ~GroupSolveIteration() = default;

  void Iterate(system::System &system) override;

  bool is_storing_angular_solution() const {
    return is_storing_angular_solution_;
  }

  EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map() const {
    return angular_solution_ptr_map_;
  }

  GroupSolver* group_solver_ptr() const {
    return group_solver_ptr_.get();
  }

  ConvergenceChecker* convergence_checker_ptr() const {
    return convergence_checker_ptr_.get();
  }

  MomentCalculator* moment_calculator_ptr() const {
    return moment_calculator_ptr_.get();
  }

  MomentMapConvergenceChecker* moment_map_convergence_checker_ptr() const {
    return moment_map_convergence_checker_ptr_.get();
  }

  std::shared_ptr<GroupSolution> group_solution_ptr() const {
    return group_solution_ptr_;
  }

  Reporter* reporter_ptr() const {
    return reporter_ptr_.get();
  }

 protected:
  virtual void PerformPerGroup(system::System& system, const int group);
  virtual void SolveGroup(const int group, system::System &system);
  virtual void StoreAngularSolution(system::System& system, const int group);
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
  std::shared_ptr<Reporter> reporter_ptr_ = nullptr;
  std::unique_ptr<MomentMapConvergenceChecker>
      moment_map_convergence_checker_ptr_ = nullptr;
  bool is_storing_angular_solution_ = false;
  EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map_;
};

} // namespace group

} // namespace iteration



} //namespace bart

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
