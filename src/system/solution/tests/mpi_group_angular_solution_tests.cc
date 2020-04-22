#include "system/solution/mpi_group_angular_solution.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

using ::testing::Ref, ::testing::Eq;

class SolutionMPIAngularTests : public ::testing::Test {
 protected:
  SolutionMPIAngularTests() : test_solution(test_angles_) {};

  // Object to be tested
  system::solution::MPIGroupAngularSolution test_solution;

  // Test parameters
  static constexpr int test_angles_ = 4;
};

TEST_F(SolutionMPIAngularTests, Constructor) {
  EXPECT_EQ(test_solution.total_angles(), test_angles_);
  EXPECT_EQ(test_solution.solutions().size(), test_angles_);

  for (int angle = 0; angle < test_angles_; ++angle) {
    EXPECT_NO_THROW(test_solution.solutions().at(angle));
  }

}

TEST_F(SolutionMPIAngularTests, OperatorBraketsPair) {

  const auto& const_test_solution = test_solution;

  for (int angle = 0; angle < test_angles_; ++angle) {
      EXPECT_THAT(test_solution[angle],
                  Ref(test_solution.solutions().at(angle)));
      EXPECT_THAT(const_test_solution[angle],
                  Ref(const_test_solution.solutions().at(angle)));
      EXPECT_THAT(test_solution.GetSolution(angle),
                  Ref(test_solution.solutions().at(angle)));
  }
}

} // namespace