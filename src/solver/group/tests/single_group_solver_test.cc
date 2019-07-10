#include "solver/group/single_group_solver.h"

#include <memory>

#include "solver/tests/linear_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class SolverGroupSingleGroupSolverTest : public ::testing::Test {
 protected:

  using LinearSolver = solver::LinearMock;

  // Supporting mocks
  std::unique_ptr<LinearSolver> linear_solver_ptr_;

  LinearSolver* linear_solver_obs_ptr_;

  void SetUp() override;
};

void SolverGroupSingleGroupSolverTest::SetUp() {
  linear_solver_ptr_ = std::make_unique<LinearSolver>();
  linear_solver_obs_ptr_ = linear_solver_ptr_.get();
}



TEST_F(SolverGroupSingleGroupSolverTest, Constructor) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

  auto test_ptr = dynamic_cast<LinearSolver*>(test_solver.linear_solver_ptr());

  EXPECT_NE(test_ptr, nullptr);
}

} // namespace
