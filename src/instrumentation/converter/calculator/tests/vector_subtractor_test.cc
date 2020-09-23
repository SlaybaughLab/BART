#include "instrumentation/converter/calculator/vector_subtractor.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace {

namespace instrumentation = bart::instrumentation;

class InstrumentationConverterCalculatorVectorSubtractorTest
 : public ::testing::Test {
 public:
  using VectorSubtractor = instrumentation::converter::calculator::VectorSubtractor;
  using DealiiVector = VectorSubtractor::DealiiVector;
  void SetUp() override;
};

void InstrumentationConverterCalculatorVectorSubtractorTest::SetUp() {

}

TEST_F(InstrumentationConverterCalculatorVectorSubtractorTest, Dummy) {
  EXPECT_TRUE(false);
}



} // namespace bart
