#include "instrumentation/converter/to_string/string_color_pair_to_string.h"

#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

using namespace bart;

class InstrumentationConverterStringColorPairToStringTest
 : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::StringColorPairToString;
  using Color = utility::Color;
  using PairType = std::pair<std::string, Color>;


  std::string default_output{"${COLOR_CODE}${STRING}${COLOR_RESET_CODE}"};
  std::string GetExpectedOutput(const PairType, const ConverterType&,
                                const std::string) const;
};

std::string
InstrumentationConverterStringColorPairToStringTest::GetExpectedOutput(
    const PairType to_convert, const ConverterType& converter,
    const std::string output_format) const {
  using OutputTerm = ConverterType::OutputTerm;
  const auto output_term_to_string_map = converter.output_term_to_string_map();
  std::string output = output_format;

  std::string string_string = output_term_to_string_map.at(OutputTerm::kString);
  if (auto index = output.find(string_string); index != std::string::npos)
    output.replace(index, string_string.size(), to_convert.first);
  std::string color_string = output_term_to_string_map.at(
      OutputTerm::kColorCode);
  if (auto index = output.find(color_string); index != std::string::npos) {
    output.replace(index, color_string.size(),
                   utility::to_string(to_convert.second));
  }
  std::string reset_string = output_term_to_string_map.at(
      OutputTerm::kColorReset);
  if (auto index = output.find(reset_string); index != std::string::npos) {
    output.replace(index, reset_string.size(),
                   utility::to_string(Color::kReset));
  }
  return output;
}

TEST_F(InstrumentationConverterStringColorPairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_ptr_;
  EXPECT_NO_THROW({test_ptr_ = std::make_unique<ConverterType>();});
  EXPECT_EQ(test_ptr_->output_format(), default_output);
}

TEST_F(InstrumentationConverterStringColorPairToStringTest, ConvertDefault) {
  ConverterType test_converter;
  PairType test_pair{"Color this string red", Color::kRed};
  EXPECT_EQ(test_converter.Convert(test_pair),
            GetExpectedOutput(test_pair, test_converter, default_output));
}

TEST_F(InstrumentationConverterStringColorPairToStringTest,
       ConvertNewOutputFormat) {
  using OutputTerm = ConverterType::OutputTerm;
  ConverterType test_converter;
  PairType test_pair{"Color this string red", Color::kRed};
  auto new_format = test_converter.SetOutputFormat(
      {"New format, colored string: ", OutputTerm::kColorCode,
       OutputTerm::kString,
       OutputTerm::kColorReset, "\n"});
  EXPECT_EQ(test_converter.Convert(test_pair),
            GetExpectedOutput(test_pair, test_converter, new_format));
}


} // namespace
