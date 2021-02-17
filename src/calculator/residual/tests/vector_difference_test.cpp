#include "calculator/residual/vector_difference.hpp"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.hpp"
#include "test_helpers/test_helper_functions.h"

#include "calculator/residual/tests/vector_difference_mock.hpp" // mock to verify good header

namespace  {

using namespace bart;

class CalculatorResidualVectorDifferenceTest : public ::testing::Test {
 public:
  using ResidualCalculator = calculator::residual::VectorDifference;
  using Vector = dealii::Vector<double>;

  // Test object
  ResidualCalculator test_calculator_;

  // Test supporting objects
  Vector minuend_, subtrahend_, expected_residual_, expected_weighed_residual_;

  // Test parameters
  const int vector_size_{ test_helpers::RandomInt(10, 20) };
  const double weight_{ test_helpers::RandomDouble(-100, 100) };

  auto SetUp() -> void override;

};

auto CalculatorResidualVectorDifferenceTest::SetUp() -> void {
  minuend_.reinit(vector_size_);
  subtrahend_.reinit(vector_size_);
  expected_residual_.reinit(vector_size_);
  expected_weighed_residual_.reinit(vector_size_);

  for (int i = 0; i < vector_size_; ++i) {
    const double minuend_value{ test_helpers::RandomDouble(-100, 100) };
    const double subtrahend_value{ test_helpers::RandomDouble(-100, 100) };

    minuend_[i] = minuend_value;
    subtrahend_[i] = subtrahend_value;
    expected_residual_[i] = minuend_value - subtrahend_value;
    expected_weighed_residual_[i] = weight_*(minuend_value - subtrahend_value);
  }
}

// Should calculate the correct residual with no weight
TEST_F(CalculatorResidualVectorDifferenceTest, CalculateResidualNoWeight) {
  const auto calculated_residual = test_calculator_.CalculateResidual(minuend_, subtrahend_);
  ASSERT_EQ(calculated_residual.size(), this->vector_size_);
  EXPECT_TRUE(test_helpers::AreEqual(calculated_residual, this->expected_residual_));
}

// Should calculate the correct weighted residual
TEST_F(CalculatorResidualVectorDifferenceTest, CalculateResidualWeighed) {
  const auto calculated_residual = test_calculator_.CalculateResidual(minuend_, subtrahend_, weight_);
  ASSERT_EQ(calculated_residual.size(), this->vector_size_);
  EXPECT_TRUE(test_helpers::AreEqual(calculated_residual, expected_weighed_residual_));
}

} // namespace
