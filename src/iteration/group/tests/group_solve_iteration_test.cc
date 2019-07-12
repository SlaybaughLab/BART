#include "iteration/group/group_solve_iteration.h"

#include <memory>

#include "solver/group/tests/single_group_solver_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IterationGroupSolveIterationTest : public ::testing::Test {
 public:
  using GroupSolver = solver::group::SingleGroupSolverMock;

  // Mock objects
  std::unique_ptr<GroupSolver> single_group_solver_ptr_;

  // Observing pointers
  GroupSolver* single_group_obs_ptr_ = nullptr;

  void SetUp() override;
};

void IterationGroupSolveIterationTest::SetUp() {
  single_group_solver_ptr_ = std::make_unique<GroupSolver>();
  single_group_obs_ptr_ = single_group_solver_ptr_.get();
}

TEST_F(IterationGroupSolveIterationTest, Constructor) {
  iteration::group::GroupSolveIteration test_iterator(
      std::move(single_group_solver_ptr_));

  auto single_group_test_ptr = dynamic_cast<GroupSolver*>(
      test_iterator.group_solver());

  EXPECT_NE(nullptr, single_group_test_ptr);
}

} // namespace