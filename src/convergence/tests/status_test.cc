#include "convergence/status.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace convergence = bart::convergence;
namespace test_helpers = bart::test_helpers;

class ConvergenceStatusOperatorsTest : public ::testing::Test {
 public:
  convergence::Status test_status_one, test_status_two;
  void SetUp() override;
  const int max_iterations{100};
};

void ConvergenceStatusOperatorsTest::SetUp() {
  test_status_one.max_iterations = max_iterations;
  test_status_one.iteration_number = test_helpers::RandomInt(0, max_iterations);
  test_status_one.is_complete = false;
  test_status_one.failed_index = std::nullopt;
  test_status_one.delta = test_helpers::RandomDouble(1e-16, 1e-4);

  test_status_two.max_iterations =   test_status_one.max_iterations;
  test_status_two.iteration_number = test_status_one.iteration_number;
  test_status_two.is_complete =      test_status_one.is_complete;
  test_status_two.failed_index =     test_status_one.failed_index;
  test_status_two.delta =            test_status_one.delta;
}

TEST_F(ConvergenceStatusOperatorsTest, CopyAssignment) {
  convergence::Status new_status = test_status_one;
  EXPECT_TRUE(new_status == test_status_one);
}

TEST_F(ConvergenceStatusOperatorsTest, CopyConstructor) {
  convergence::Status new_status(test_status_one);
  EXPECT_TRUE(new_status == test_status_one);
}

TEST_F(ConvergenceStatusOperatorsTest, Equivalence) {
  EXPECT_TRUE(test_status_one == test_status_two);
  EXPECT_FALSE(test_status_one != test_status_two);
}

TEST_F(ConvergenceStatusOperatorsTest, NonEquivalenceIterationNumber) {
  ++test_status_two.iteration_number;
  EXPECT_FALSE(test_status_one == test_status_two);
  EXPECT_TRUE(test_status_one != test_status_two);
}

TEST_F(ConvergenceStatusOperatorsTest, NonEquivalenceMaxIterations) {
  ++test_status_two.max_iterations;
  EXPECT_FALSE(test_status_one == test_status_two);
  EXPECT_TRUE(test_status_one != test_status_two);
}

TEST_F(ConvergenceStatusOperatorsTest, NonEquivalenceIsComplete) {
  test_status_two.is_complete = !test_status_one.is_complete;
  EXPECT_FALSE(test_status_one == test_status_two);
  EXPECT_TRUE(test_status_one != test_status_two);
}

TEST_F(ConvergenceStatusOperatorsTest, NonEquivalenceFailedIndex) {
  test_status_two.failed_index = 3;
  EXPECT_FALSE(test_status_one == test_status_two);
  EXPECT_TRUE(test_status_one != test_status_two);
}

TEST_F(ConvergenceStatusOperatorsTest, NonEquivalenceDelta) {
  test_status_two.delta = std::nullopt;
  EXPECT_FALSE(test_status_one == test_status_two);
  EXPECT_TRUE(test_status_one != test_status_two);
}

} // namespace