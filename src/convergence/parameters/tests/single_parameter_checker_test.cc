#include "convergence/parameters/single_parameter_checker.h"

#include <vector>

#include "convergence/tests/single_checker_test.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

using namespace bart;

class ConvergenceSingleParameterCheckerTest :
 public bart::convergence::testing::SingleCheckerTest<double>{

 protected:
  bart::convergence::parameters::SingleParameterChecker test_checker{1e-6};
  std::vector<double> test_values;
  void SetUp() override;
};
void ConvergenceSingleParameterCheckerTest::SetUp() {
  // Get four random doubles between 0.01 and 10
  test_values = test_helpers::RandomVector(4, 0.01, 10);
}

TEST_F(ConvergenceSingleParameterCheckerTest, BaseClassTests) {
  TestBaseMethods(&test_checker);
}

TEST_F(ConvergenceSingleParameterCheckerTest, SameValue) {
  for (auto const val : test_values) {
    EXPECT_TRUE(test_checker.CheckIfConverged(val, val));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val));
  }
}

TEST_F(ConvergenceSingleParameterCheckerTest, OneThreshold) {
  for (auto val : test_values) {
    auto val_one_more = val + 0.99*test_checker.max_delta()*val;
    auto val_one_less = val - 0.99*test_checker.max_delta()*val;

    EXPECT_TRUE(test_checker.CheckIfConverged(val, val_one_more));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val_one_more));
    EXPECT_TRUE(test_checker.CheckIfConverged(val, val_one_less));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val_one_less));
  }
}

TEST_F(ConvergenceSingleParameterCheckerTest, MoreThanOneThreshold) {

  std::vector<double> multiples{1.1, 1.5, 2, 10};

  for (auto val : test_values) {
    for (auto const multiple : multiples){
      auto val_more = val + multiple*test_checker.max_delta()*val;
      auto val_less = val - multiple*test_checker.max_delta()*val;

      EXPECT_FALSE(test_checker.CheckIfConverged(val, val_more));
      EXPECT_FALSE(test_checker.CheckIfConverged(-val, -val_more));
      EXPECT_FALSE(test_checker.CheckIfConverged(val, val_less));
      EXPECT_FALSE(test_checker.CheckIfConverged(-val, -val_less));
    }
  }
}



} // namespace

