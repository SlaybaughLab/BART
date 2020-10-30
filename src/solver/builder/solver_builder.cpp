#include "solver/builder/solver_builder.hpp"

#include "solver/group/factory.hpp"
#include "solver/linear/factory.hpp"

namespace bart::solver::builder {

template <>
auto SolverBuilder::BuildSolver(const SolverName name, const int max_iterations, const double convergence_tolerance)
-> std::unique_ptr<group::SingleGroupSolverI> {
  // Build linear solver
  std::unique_ptr<linear::LinearI> linear_solver_ptr;
  switch (name) {
    case SolverName::kDefaultGMRESGroupSolver: {
      linear_solver_ptr = std::move(linear::LinearIFactory<int, double>::get()
                                        .GetConstructor(linear::LinearSolverName::kGMRES)
                                            (max_iterations, convergence_tolerance));
    }
  }

  // Build group solver
  std::unique_ptr<group::SingleGroupSolverI> return_ptr;
  switch (name) {
    case SolverName::kDefaultGMRESGroupSolver: {
      return solver::group::SingleGroupSolverIFactory<std::unique_ptr<linear::LinearI>>::get()
          .GetConstructor(solver::group::GroupSolverName::kDefaultImplementation)
              (std::move(linear_solver_ptr));
    }
  }
  return nullptr;
}

template <>
auto SolverBuilder::BuildSolver(const SolverName name) -> std::unique_ptr<group::SingleGroupSolverI> {
  switch (name) {
    case SolverName::kDefaultGMRESGroupSolver: {
      return BuildSolver(SolverName::kDefaultGMRESGroupSolver, 100, 1e-10);
    }
  }
  return nullptr;
}

} // namespace bart::solver::builder
