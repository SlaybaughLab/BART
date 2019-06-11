#include "system/solution/mpi_angular.h"

#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

class SolutionMPIAngularTests : public ::testing::Test {
 protected:
  SolutionMPIAngularTests() : test_solution(test_groups_, test_angles_) {};

  // Object to be tested
  system::solution::MPIAngular test_solution;

  // Test parameters
  static constexpr int test_groups_ = 2;
  static constexpr int test_angles_ = 4;
};

TEST_F(SolutionMPIAngularTests, Constructor) {
  EXPECT_EQ(test_solution.total_angles(), test_angles_);
  EXPECT_EQ(test_solution.total_groups(), test_groups_);
  EXPECT_EQ(test_solution.solutions().size(), test_groups_*test_angles_);

  for (int group = 0; group < test_groups_; ++group) {
    for (int angle = 0; angle < test_angles_; ++angle) {
      EXPECT_NO_THROW(test_solution.solutions().at({group, angle}));
    }
  }
}

} // namespace