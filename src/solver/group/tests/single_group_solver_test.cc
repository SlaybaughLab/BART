#include "solver/group/single_group_solver.h"

#include <memory>

#include "system/system.h"
#include "system/solution/tests/mpi_group_angular_solution_mock.h"
#include "system/terms/tests/linear_term_mock.h"
#include "system/terms/tests/bilinear_term_mock.h"
#include "solver/linear/tests/linear_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return;
using ::testing::ReturnRef, ::testing::_;
using ::testing::Pointee, ::testing::Ref;

class SolverGroupSingleGroupSolverTest :
    public ::testing::Test,
    public bart::testing::DealiiTestDomain<2> {
 protected:

  using LinearSolver = solver::linear::LinearMock;
  using LeftHandSide = system::terms::BilinearTermMock;
  using RightHandSide = system::terms::LinearTermMock;
  using GroupSolution = NiceMock<system::solution::MPIGroupAngularSolutionMock>;

  // SUpporting objects
  system::System test_system_;
  GroupSolution solution_;

  // Supporting mocks
  std::unique_ptr<LinearSolver> linear_solver_ptr_;

  // Mock Observing pointers
  LinearSolver* linear_solver_obs_ptr_;
  RightHandSide* rhs_obs_ptr_;
  LeftHandSide* lhs_obs_ptr_;

  // test parameters
  const int total_angles_ = 2;
  const int test_group_ = 2;

  void SetUp() override;
};

void SolverGroupSingleGroupSolverTest::SetUp() {
  SetUpDealii();
  linear_solver_ptr_ = std::make_unique<LinearSolver>();
  linear_solver_obs_ptr_ = linear_solver_ptr_.get();

  auto rhs_ptr_ = std::make_unique<RightHandSide>();
  auto lhs_ptr_ = std::make_unique<LeftHandSide>();

  rhs_obs_ptr_ = rhs_ptr_.get();
  lhs_obs_ptr_ = lhs_ptr_.get();

  test_system_.right_hand_side_ptr_ = std::move(rhs_ptr_);
  test_system_.left_hand_side_ptr_ = std::move(lhs_ptr_);

  ON_CALL(solution_, total_angles())
      .WillByDefault(Return(total_angles_));
}

TEST_F(SolverGroupSingleGroupSolverTest, Constructor) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

  auto test_ptr = dynamic_cast<LinearSolver*>(test_solver.linear_solver_ptr());

  EXPECT_NE(test_ptr, nullptr);
}

TEST_F(SolverGroupSingleGroupSolverTest, SolveGroupOperation) {
  solver::group::SingleGroupSolver test_solver(std::move(linear_solver_ptr_));

  /* Test objects. These are objects that will be returned by our mock functions
   * and used to verify correct mediation of all the classes. In most cases they
   * do not need to include any actual information, we will just verify they are
   * correctly passed by their reference.
   */
  std::vector<system::MPIVector> solution_vectors_(total_angles_);
  std::vector<std::shared_ptr<system::MPISparseMatrix>> lhs_matrices_(total_angles_);
  std::vector<std::shared_ptr<system::MPIVector>> rhs_vectors_(total_angles_);

  // Expect to retrieve total angles (to access all angles for the test group)
  EXPECT_CALL(solution_, total_angles())
      .WillOnce(Return(total_angles_));

  // Expectations for each angle:
  for (int angle = 0; angle < total_angles_; ++angle) {
    system::Index index{test_group_, angle};

    rhs_vectors_[angle] = std::make_shared<system::MPIVector>();
    lhs_matrices_[angle] = std::make_shared<system::MPISparseMatrix>();
    lhs_matrices_[angle]->reinit(matrix_1);
    lhs_matrices_[angle]->copy_from(matrix_1);

    // Expect to retrieve each solution
    EXPECT_CALL(solution_, BracketOp(angle))
        .WillOnce(ReturnRef(solution_vectors_[angle]));
    // Expect to retrieve each LHS
    EXPECT_CALL(*lhs_obs_ptr_, GetFullTermPtr(index))
        .WillOnce(Return(lhs_matrices_[angle]));
    // Expect to retrieve each RHS
    EXPECT_CALL(*rhs_obs_ptr_, GetFullTermPtr(index))
        .WillOnce(Return(rhs_vectors_[angle]));

    EXPECT_CALL(*linear_solver_obs_ptr_, Solve(
        lhs_matrices_[angle].get(),
        Pointee(solution_vectors_[angle]),
        rhs_vectors_[angle].get(),
        _));
  }

  test_solver.SolveGroup(test_group_, test_system_, solution_);
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


} // namespace
