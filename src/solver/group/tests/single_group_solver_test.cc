#include "solver/group/single_group_solver.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class SolverGroupSingleGroupSolverTest : public ::testing::Test {
 protected:
};

TEST_F(SolverGroupSingleGroupSolverTest, Constructor) {
  solver::group::SingleGroupSolver test_solver;
}

} // namespace
