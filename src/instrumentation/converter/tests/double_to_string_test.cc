#include "instrumentation/converter/double_to_string.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterDoubleToStringTest : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::DoubleToString;
  const double test_value_{test_helpers::RandomDouble(0, 100)};
  std::string GetExpectedOutput(const double, ConverterType&, std::string) const;
};

std::string InstrumentationConverterDoubleToStringTest::GetExpectedOutput(
    const double to_convert, ConverterType& converter,
    std::string overload_format = "") const {
  using OutputTerm = ConverterType::OutputTerm;
}

TEST_F(InstrumentationConverterDoubleToStringTest, Constructor) {
  EXPECT_NO_THROW({
    ConverterType test_converter;
  });
}

} // namespace
