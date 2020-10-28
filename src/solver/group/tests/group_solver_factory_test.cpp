#include "solver/group/factory.hpp"

#include "solver/group/single_group_solver.h"
#include "solver/linear/tests/linear_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace solver = bart::solver;

TEST(SolverGroupFactoryTests, SingleGroupSolverDefaultImplementation) {
  using SolverName = solver::group::GroupSolverName;
  using ExpectedType = solver::group::SingleGroupSolver;
  using LinearSolver = solver::LinearI;

  auto group_solver_ptr =
      solver::group::SingleGroupSolverIFactory<std::unique_ptr<LinearSolver>>::get()
          .GetConstructor(SolverName::kDefaultImplementation)(
              std::make_unique<solver::LinearMock>());
  ASSERT_NE(group_solver_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(group_solver_ptr.get()), nullptr);
}

} // namespace
