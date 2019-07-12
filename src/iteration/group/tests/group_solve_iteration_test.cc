#include "iteration/group/group_solve_iteration.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class IterationGroupSolveIterationTest : public ::testing::Test {

};

TEST_F(IterationGroupSolveIterationTest, Constructor) {
  iteration::group::GroupSolveIteration group_iterator;
}

} // namespace