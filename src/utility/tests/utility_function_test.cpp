#include "../utility_functions.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"

class TestUtilityFunctions : public ::testing::Test {
 protected:
};

// TEST_F(TestUtilityFunctions, VectorTest) {
//   std::vector<double> numbers =
//       {0.01, 0.001, 0.0001, 0.000001, 0.00000000001};
//   double sum = bart::utility::PreciseSum(numbers);
//   EXPECT_FLOAT_EQ(sum, 0.011101);
// }

TEST_F(TestUtilityFunctions, VectorTest) {
  std::vector<double> numbers =
      {10e100, 1.0, -10e100, 1.0};
  double sum = bart::utility::PreciseSum(numbers);
  EXPECT_FLOAT_EQ(sum, 2);
}

TEST_F(TestUtilityFunctions, SecondVectorTest) {
  std::vector<double> numbers =
      {1.0, 10e100, 1.0, -10e100};
  double sum = bart::utility::PreciseSum(numbers);
  EXPECT_FLOAT_EQ(sum, 2);
}

TEST_F(TestUtilityFunctions, ArrayTest) {
  std::array<double, 5> numbers =
      {0.01, 0.001, 0.0001, 0.000001, 0.00000000001};
  double sum = bart::utility::PreciseSum(numbers);
  EXPECT_FLOAT_EQ(sum, 0.011101);
}




