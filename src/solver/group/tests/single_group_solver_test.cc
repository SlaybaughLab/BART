#include "solver/group/single_group_solver.h"

#include <memory>

#include "system/system.h"
#include "system/solution/tests/mpi_angular_mock.h"
#include "solver/tests/linear_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;

class SolverGroupSingleGroupSolverTest : public ::testing::Test {
 protected:

  using LinearSolver = solver::LinearMock;
  using GroupSolution = NiceMock<system::solution::MPIAngularMock>;

  // SUpporting objects
  system::System test_system_;
  GroupSolution solution_;

  // Supporting mocks
  std::unique_ptr<LinearSolver> linear_solver_ptr_;


  LinearSolver* linear_solver_obs_ptr_;

  // test parameters
  const int total_angles_ = 2;
  const int total_groups_ = 4;

  void SetUp() override;
};

void SolverGroupSingleGroupSolverTest::SetUp() {
  linear_solver_ptr_ = std::make_unique<LinearSolver>();
  linear_solver_obs_ptr_ = linear_solver_ptr_.get();

  ON_CALL(solution_, total_angles())
      .WillByDefault(Return(total_angles_));
  ON_CALL(solution_, total_groups())
      .WillByDefault(Return(total_groups_));
}

TEST_F(SolverGroupSingleGroupSolverTest, Constructor) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

  auto test_ptr = dynamic_cast<LinearSolver*>(test_solver.linear_solver_ptr());

  EXPECT_NE(test_ptr, nullptr);
}

TEST_F(SolverGroupSingleGroupSolverTest, SolveGroupBadAngles) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

  std::array<int, 2> bad_angles = {-1, 0};
  for (const int angle : bad_angles) {
    EXPECT_CALL(solution_, total_angles())
        .WillOnce(Return(angle));
    EXPECT_ANY_THROW(test_solver.SolveGroup(0, test_system_, solution_));
  }
}

TEST_F(SolverGroupSingleGroupSolverTest, SolveGroupBadGroup) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));
  EXPECT_ANY_THROW(test_solver.SolveGroup(4, test_system_, solution_));
  EXPECT_ANY_THROW(test_solver.SolveGroup(-1, test_system_, solution_));
}

TEST_F(SolverGroupSingleGroupSolverTest, SolveGroup) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

//  EXPECT_CALL(solution_, total_angles())
//      .WillOnce(Return(total_angles_));
//

}

} // namespace
