#include "instrumentation/converter/double_to_string.h"

#include <iomanip>
#include <sstream>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterDoubleToStringTest : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::DoubleToString;
  const double test_value_{test_helpers::RandomDouble(0, 100)};
  std::string GetExpectedOutput(const double, ConverterType&, std::string,
                                const int) const;
};

std::string InstrumentationConverterDoubleToStringTest::GetExpectedOutput(
    const double to_convert, ConverterType& converter,
    std::string overload_format = "", const int precision = -1) const {
  using OutputTerm = ConverterType::OutputTerm;
  const auto output_term_to_string_map = converter.output_term_to_string_map();
  int output_precision = precision;
  std::string output = overload_format;
  if (overload_format == "")
    output = converter.output_format();
  if (precision == -1)
    output_precision = converter.precision();

  std::string value_string = output_term_to_string_map.at(OutputTerm::kValue);
  auto value_index = output.find(value_string);
  if (value_index != std::string::npos) {
    std::ostringstream string_stream;
    string_stream << std::fixed
                  << std::setprecision(output_precision)
                  << to_convert;
    output.replace(value_index, value_string.size(),
                   string_stream.str());
  }

  return output;
}

TEST_F(InstrumentationConverterDoubleToStringTest, Constructor) {
  EXPECT_NO_THROW({
    ConverterType test_converter;
    EXPECT_EQ(test_converter.precision(), 2);
    EXPECT_EQ(test_converter.output_format(), "${VALUE}\n");
  });
}

TEST_F(InstrumentationConverterDoubleToStringTest, DefaultStringAndPrecision) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.Convert(test_value_),
            GetExpectedOutput(test_value_, test_converter));
}

TEST_F(InstrumentationConverterDoubleToStringTest, SetPrecision) {
  ConverterType test_converter;
  const int new_precision = test_helpers::RandomDouble(3, 10);
  EXPECT_EQ(test_converter.set_precision(new_precision).Convert(test_value_),
            GetExpectedOutput(test_value_, test_converter, "", new_precision));
}

TEST_F(InstrumentationConverterDoubleToStringTest, ChangeString) {
  ConverterType test_converter;
  using OutputTerm = ConverterType::OutputTerm;
  const auto value_string = test_converter.output_term_to_string_map()
      .at(OutputTerm::kValue);

  std::string new_format = "Value is: " + value_string + "\n";
  auto returned_format = test_converter.SetOutputFormat({
    "Value is: ", OutputTerm::kValue, "\n"});
  EXPECT_EQ(returned_format, new_format);
  EXPECT_EQ(test_converter.Convert(test_value_),
            GetExpectedOutput(test_value_, test_converter, new_format));
}

} // namespace
