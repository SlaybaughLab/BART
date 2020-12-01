#include "instrumentation/converter/calculator/vector_subtractor.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

namespace instrumentation = bart::instrumentation;
namespace test_helpers = bart::test_helpers;

class InstrumentationConverterCalculatorVectorSubtractorTest
 : public ::testing::TestWithParam<bool> {
 public:
  using VectorSubtractor = instrumentation::converter::calculator::VectorSubtractor;
  using DealiiVector = VectorSubtractor::DealiiVector;
  using AbsoluteValue = instrumentation::converter::calculator::AbsoluteValue;

  std::unique_ptr<VectorSubtractor> test_subtractor_ptr_;
  const int vector_length_{test_helpers::RandomInt(10, 20)};
  DealiiVector minuend_, subtrahend_, difference_; // difference = minuend - subtrahend
  const bool calculate_absolute_value_{GetParam()};
  void SetUp() override;
};

INSTANTIATE_TEST_SUITE_P(AbsoluteValueSubtraction,
                         InstrumentationConverterCalculatorVectorSubtractorTest,
                         ::testing::Bool());

void InstrumentationConverterCalculatorVectorSubtractorTest::SetUp() {
  minuend_.reinit(vector_length_);
  subtrahend_.reinit(vector_length_);
  difference_.reinit(vector_length_);

  for (int i = 0; i < vector_length_; ++i) {
    double minuend_value{test_helpers::RandomDouble(-100, 100)};
    double subtrahend_value{test_helpers::RandomDouble(-100, 100)};

    minuend_[i] = minuend_value;
    subtrahend_[i] = subtrahend_value;
    double difference_value = minuend_value - subtrahend_value;

    if (calculate_absolute_value_) {
      difference_value = std::abs(difference_value);
    }

    difference_[i] = difference_value;
  }
  test_subtractor_ptr_ = std::make_unique<VectorSubtractor>(
      minuend_,
      AbsoluteValue(calculate_absolute_value_));
}

TEST_P(InstrumentationConverterCalculatorVectorSubtractorTest, Getters) {
  EXPECT_EQ(test_subtractor_ptr_->calculate_absolute_value(),
            calculate_absolute_value_);
  ASSERT_EQ(test_subtractor_ptr_->minuend().size(), minuend_.size());
  EXPECT_EQ(test_subtractor_ptr_->minuend(), minuend_);
}

TEST_P(InstrumentationConverterCalculatorVectorSubtractorTest, ConvertBadSizes) {

  for (const int bad_size : {vector_length_ - 1, vector_length_ + 1}) {
    DealiiVector bad_subtrahend(bad_size);
    for (dealii::Vector<double>::size_type i = 0; i < bad_subtrahend.size(); ++i)
      bad_subtrahend[i] = test_helpers::RandomDouble(-100, 100);
    EXPECT_ANY_THROW(test_subtractor_ptr_->Convert(bad_subtrahend))
              << "Failed to throw of bad vector size: " << bad_size;
  }
}

TEST_P(InstrumentationConverterCalculatorVectorSubtractorTest, Convert) {
  auto result = test_subtractor_ptr_->Convert(subtrahend_);
  ASSERT_EQ(result.size(), difference_.size());
  EXPECT_EQ(result, difference_);
}






} // namespace bart
