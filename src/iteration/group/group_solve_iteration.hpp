#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_HPP_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_HPP_

#include <memory>

#include "convergence/iteration_completion_checker_i.hpp"
#include "instrumentation/port.hpp"
#include "iteration/group/group_solve_iteration_i.hpp"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"
#include "solver/group/single_group_solver_i.h"
#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/solution/solution_types.h"
#include "utility/has_dependencies.h"

namespace bart::iteration::group {

//! Data ports for group iterations.
namespace data_ports {
struct GroupConvergenceStatus;
struct Status;
struct NumberOfIterations;
//! Data port for the status of convergence.
using ConvergenceStatusPort = instrumentation::Port<convergence::Status, GroupConvergenceStatus>;
using NumberOfIterationsPort = instrumentation::Port<double, NumberOfIterations>;
//! Data port for general strings.
using StatusPort = instrumentation::Port<std::string, Status>;
}

/*! \brief Default implementation for group solve iterations.
 *
 * This base class provides the basis process for converging all groups. The Iterate method does the following:
 *
 * 1. Updates system previous moments.
 * 2. For each group:
 *   a. Saves the current group scalar flux.
 *   b. Updates the system for the current group.
 *   c. Solves the group.
 *   d. Checks for scalar flux convergence. If not converged, returns to 2.a.
 * 3. Checks that all group scalar fluxes have converged.
 *
 * The only portion that is not specified by this base class is 2b, updating the system.
 *
 */
template <int dim>
class GroupSolveIteration : public GroupSolveIterationI, public utility::HasDependencies,
                            public data_ports::ConvergenceStatusPort, public data_ports::StatusPort,
                            public data_ports::NumberOfIterationsPort {
 public:
  using GroupSolver = solver::group::SingleGroupSolverI;
  using ConvergenceChecker = convergence::IterationCompletionCheckerI<system::moments::MomentVector>;
  using MomentMapConvergenceChecker = convergence::IterationCompletionCheckerI<system::moments::MomentsMap>;
  using MomentCalculator = quadrature::calculators::SphericalHarmonicMomentsI;
  using MomentVector = system::moments::MomentVector;
  using GroupSolution = system::solution::MPIGroupAngularSolutionI;
  using EnergyGroupToAngularSolutionPtrMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  using Subroutine = iteration::subroutine::SubroutineI;
  using System = system::System;

  // Data ports
  using data_ports::ConvergenceStatusPort::Expose, data_ports::ConvergenceStatusPort::AddInstrument;
  using data_ports::StatusPort::Expose, data_ports::StatusPort::AddInstrument;

  GroupSolveIteration(std::unique_ptr<GroupSolver>, std::unique_ptr<ConvergenceChecker>,
      std::unique_ptr<MomentCalculator>, std::shared_ptr<GroupSolution>,
      std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr = nullptr);

  void Iterate(System &system) override;
  auto UpdateThisAngularSolutionMap(EnergyGroupToAngularSolutionPtrMap to_update) -> GroupSolveIteration& override {
    is_storing_angular_solution_ = true;
    angular_solution_ptr_map_ = to_update;
    return *this;
  }

  [[nodiscard]] auto is_storing_angular_solution() const -> bool { return is_storing_angular_solution_; }
  [[nodiscard]] auto angular_solution_ptr_map() const { return angular_solution_ptr_map_; }

  auto AddPostIterationSubroutine(std::unique_ptr<Subroutine> subroutine_ptr) -> GroupSolveIteration<dim>& override {
    post_iteration_subroutine_ptr_ = std::move(subroutine_ptr);
    return *this; };

  auto group_solver_ptr() const { return group_solver_ptr_.get(); }
  auto convergence_checker_ptr() const { return convergence_checker_ptr_.get(); }
  auto moment_calculator_ptr() const { return moment_calculator_ptr_.get(); }
  auto moment_map_convergence_checker_ptr() const { return moment_map_convergence_checker_ptr_.get(); }
  [[nodiscard]] auto group_solution_ptr() const -> std::shared_ptr<GroupSolution> { return group_solution_ptr_; }
  auto post_iteration_subroutine_ptr() const { return post_iteration_subroutine_ptr_.get(); }
 protected:
  virtual auto PerformPerGroup(System& system, int group) -> void;
  virtual auto SolveGroup(int group, System &system) -> void;
  virtual auto StoreAngularSolution(System& system, int group) -> void;
  virtual auto GetScalarFlux(int group, System& system) -> MomentVector;
  virtual auto CheckConvergence(const MomentVector& current_iteration,
                                const MomentVector& previous_iteration) -> convergence::Status;
  virtual auto UpdateSystem(System& system, int group, int angle) -> void = 0;
  virtual auto UpdateCurrentMoments(System &system, int group) -> void;
  virtual auto ExposeIterationData(system::System&) -> void {};

  std::unique_ptr<GroupSolver> group_solver_ptr_{ nullptr };
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_{ nullptr };
  std::unique_ptr<MomentCalculator> moment_calculator_ptr_{ nullptr };
  std::shared_ptr<GroupSolution> group_solution_ptr_{ nullptr };
  std::unique_ptr<MomentMapConvergenceChecker> moment_map_convergence_checker_ptr_{ nullptr };
  std::unique_ptr<Subroutine> post_iteration_subroutine_ptr_{ nullptr };
  bool is_storing_angular_solution_{ false };
  EnergyGroupToAngularSolutionPtrMap angular_solution_ptr_map_{};
};

} // namespace bart::iteration::group

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_HPP_
