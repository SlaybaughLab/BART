#include "convergence/parameters/single_parameter_checker.hpp"

#include <vector>

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

using namespace bart;

class ConvergenceSingleParameterCheckerTest : public ::testing::Test {
 protected:
  static constexpr double max_delta{ 1e-6 };
  bart::convergence::parameters::SingleParameterChecker test_checker{ max_delta };
  std::vector<double> test_values;
  auto SetUp() -> void override;
};
auto ConvergenceSingleParameterCheckerTest::SetUp() -> void {
  // Get four random doubles between 0.01 and 10
  test_values = test_helpers::RandomVector(4, 0.01, 10);
}

// Getter should return the correct value
TEST_F(ConvergenceSingleParameterCheckerTest, Getters) {
  EXPECT_EQ(test_checker.max_delta(), this->max_delta);
}

// Max delta setter should set the correct value
TEST_F(ConvergenceSingleParameterCheckerTest, MaxDeltaSetter) {
  const double max_delta_to_set{ test_helpers::RandomDouble(1e-16, 1e-8) };
  EXPECT_EQ(test_checker.max_delta(), this->max_delta);
  test_checker.SetMaxDelta(max_delta_to_set);
  EXPECT_EQ(test_checker.max_delta(), max_delta_to_set);
}

// Trying to set a negative value for the delta shoudl throw
TEST_F(ConvergenceSingleParameterCheckerTest, SetterBadValue) {
  const double max_delta_to_set{ test_helpers::RandomDouble(-1e-8, -1e-16) };
  EXPECT_EQ(test_checker.max_delta(), this->max_delta);
  EXPECT_ANY_THROW({test_checker.SetMaxDelta(max_delta_to_set);});
}

// Checking the same value should return converged
TEST_F(ConvergenceSingleParameterCheckerTest, SameValue) {
  for (auto const val : test_values) {
    EXPECT_TRUE(test_checker.CheckIfConverged(val, val));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val));
  }
}

// Being less than max delta away should return converged
TEST_F(ConvergenceSingleParameterCheckerTest, OneThreshold) {
  for (auto val : test_values) {
    const auto val_one_more = val + 0.99*test_checker.max_delta()*val;
    const auto val_one_less = val - 0.99*test_checker.max_delta()*val;

    EXPECT_TRUE(test_checker.CheckIfConverged(val, val_one_more));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val_one_more));
    EXPECT_TRUE(test_checker.CheckIfConverged(val, val_one_less));
    EXPECT_TRUE(test_checker.CheckIfConverged(-val, -val_one_less));
  }
}

// Being greater than max delta away should return false
TEST_F(ConvergenceSingleParameterCheckerTest, MoreThanOneThreshold) {
  const std::vector<double> multiples{1.1, 1.5, 2, 10};
  for (const auto val : test_values) {
    for (const auto multiple : multiples){
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

