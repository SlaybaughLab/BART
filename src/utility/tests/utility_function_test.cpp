#include "utility/utility_functions.hpp"

#include <array>
#include <vector>

#include "test_helpers/gmock_wrapper.h"

namespace  {

TEST(TestUtilityFunctions, VectorTest) {
  std::vector<double> numbers = {10e100, 1.0, -10e100, 1.0};
  EXPECT_FLOAT_EQ(bart::utility::PreciseSum(numbers), 2);
}

TEST(TestUtilityFunctions, SecondVectorTest) {
  std::vector<double> numbers = {1.0, 10e100, 1.0, -10e100};
  EXPECT_FLOAT_EQ(bart::utility::PreciseSum(numbers), 2);
}

TEST(TestUtilityFunctions, ArrayTest) {
  std::array<double, 5> numbers = {0.01, 0.001, 0.0001, 0.000001, 0.00000000001};
  EXPECT_FLOAT_EQ(bart::utility::PreciseSum(numbers), 0.011101);
}

} // namespace






