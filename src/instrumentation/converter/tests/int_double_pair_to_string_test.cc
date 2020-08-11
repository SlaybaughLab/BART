#include "instrumentation/converter/int_double_pair_to_string.h"

#include <optional>
#include <sstream>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterIntDoublePairToStringTest
    : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::IntDoublePairToString;
  const double test_double_{test_helpers::RandomDouble(20, 100)};
  const int test_int_{static_cast<int>(test_helpers::RandomDouble(0, 10))};
  std::string GetExpectedOutput(const std::pair<int, double>, const ConverterType&) const;
  std::string GetExpectedOutput(const std::pair<int, double>, const ConverterType&,
                                const std::string, const int) const;
  const int default_precision_{2};
  const std::string default_output_{"${INDEX}, ${VALUE}\n"};
};

std::string
InstrumentationConverterIntDoublePairToStringTest::GetExpectedOutput(
    const std::pair<int, double> to_convert,
    const ConverterType& converter) const {
  return GetExpectedOutput(to_convert, converter,
                           converter.output_format(),
                           converter.precision());
}
std::string
InstrumentationConverterIntDoublePairToStringTest::GetExpectedOutput(
    const std::pair<int, double> to_convert,
    const ConverterType& converter,
    const std::string format,
    const int precision) const {
  using OutputTerm = ConverterType::OutputTerm;
  const auto output_term_to_string_map = converter.output_term_to_string_map();
  std::string output = format;

  std::string index_string = output_term_to_string_map.at(OutputTerm::kIndex);
  if (std::size_t index = output.find(index_string);
      index != std::string::npos) {
    output.replace(index, index_string.size(), std::to_string(to_convert.first));
  }
  std::string value_string = output_term_to_string_map.at(OutputTerm::kValue);
  if (std::size_t index = output.find(value_string);
      index != std::string::npos) {
    std::ostringstream string_stream;
    string_stream << std::scientific << std::setprecision(precision)
                  << to_convert.second;
    output.replace(index, value_string.size(), string_stream.str());
  }
  return output;
}

TEST_F(InstrumentationConverterIntDoublePairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_converter_ptr;
  EXPECT_NO_THROW({
    test_converter_ptr = std::make_unique<ConverterType>();
  });
  EXPECT_EQ(test_converter_ptr->output_format(), default_output_);
  EXPECT_EQ(test_converter_ptr->precision(), default_precision_);
}

TEST_F(InstrumentationConverterIntDoublePairToStringTest, SetPrecision) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.precision(), default_precision_);
  const int new_precision(static_cast<int>(test_helpers::RandomDouble(3, 5)));
  test_converter.set_precision(new_precision);
  EXPECT_EQ(test_converter.precision(), new_precision);
}

TEST_F(InstrumentationConverterIntDoublePairToStringTest,
       ConvertDefaults) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.Convert({test_int_, test_double_}),
            GetExpectedOutput({test_int_, test_double_},
                              test_converter,
                              default_output_,
                              default_precision_));
}

TEST_F(InstrumentationConverterIntDoublePairToStringTest,
       ConvertNewPrecision) {
  ConverterType test_converter;
  const int new_precision{static_cast<int>(test_helpers::RandomDouble(3, 10))};
  EXPECT_EQ(test_converter
                .set_precision(new_precision)
                .Convert({test_int_, test_double_}),
            GetExpectedOutput({test_int_, test_double_},
                              test_converter,
                              default_output_,
                              new_precision));
}

TEST_F(InstrumentationConverterIntDoublePairToStringTest,
       ConvertNewOutputFormat) {
  using OutputTerm = ConverterType::OutputTerm;
  ConverterType test_converter;
  auto new_format = test_converter.SetOutputFormat(
      {"New format: index: ", OutputTerm::kIndex, " value: ", OutputTerm::kValue});
  EXPECT_EQ(test_converter
                .Convert({test_int_, test_double_}),
            GetExpectedOutput({test_int_, test_double_},
                              test_converter,
                              new_format,
                              default_precision_));
}





} // namespace
