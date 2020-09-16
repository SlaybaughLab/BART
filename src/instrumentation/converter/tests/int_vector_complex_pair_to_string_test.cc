#include "instrumentation/converter/int_vector_complex_pair_to_string.h"

#include <sstream>

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterIntVectorComplexPairToStringTest
 : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::IntVectorComplexPairToString;
  using InputType = ConverterType::InputType;
  using OutputTerm = ConverterType::OutputTerm;
  using VectorTerm = ConverterType::VectorTerm;

  std::string GetExpectedOutput(const InputType& to_convert,
                                const ConverterType& reference_converter,
                                const std::string output_format,
                                const std::string vector_format,
                                const std::string delimiter,
                                const int precision) const;

  const int default_precision_{2};
  const std::string default_output_{"${INT}, ${VECTOR}\n"};
  const std::string default_vector_entry_output_format{"(${REAL}, ${IMAG})"};
  const std::string default_delimiter_{", "};
  InputType test_input;
  void SetUp() override;
};

void InstrumentationConverterIntVectorComplexPairToStringTest::SetUp() {
  const int test_int{test_helpers::RandomInt(0, 10)};
  const int vector_size{test_helpers::RandomInt(2, 4)};
  std::vector<std::complex<double>> test_vector;
  for (int i = 0; i < vector_size; ++i) {
    test_vector.push_back({test_helpers::RandomDouble(0, 10),
                           test_helpers::RandomDouble(10, 20)});
  }
  test_input = {test_int, test_vector};
}
std::string
InstrumentationConverterIntVectorComplexPairToStringTest::GetExpectedOutput(
    const InputType& to_convert,
    const ConverterType& reference_converter,
    const std::string output_format,
    const std::string vector_format,
    const std::string delimiter,
    const int precision) const {

  const auto output_term_to_string_map =
      reference_converter.output_term_to_string_map();
  const auto vector_term_to_string_map =
      reference_converter.vector_term_to_string_map();

  std::string output = output_format;

  std::string int_string = output_term_to_string_map.at(OutputTerm::kInt);
  if (std::size_t index = output.find(int_string);
      index != std::string::npos) {
    output.replace(index, int_string.size(), std::to_string(to_convert.first));
  }
  std::string full_vector_string{""},
      index_key_string{vector_term_to_string_map.at(VectorTerm::kIndex)},
      real_key_string{vector_term_to_string_map.at(VectorTerm::kReal)},
      imag_key_string{vector_term_to_string_map.at(VectorTerm::kImag)};

  for (int i = 0; i < static_cast<int>(to_convert.second.size()); ++i) {
    std::string vector_entry{vector_format};
    if (std::size_t index = vector_entry.find(index_key_string);
        index != std::string::npos) {
      vector_entry.replace(index, index_key_string.size(), std::to_string(i));
    }
    if (std::size_t index = vector_entry.find(real_key_string);
        index != std::string::npos) {
      std::ostringstream string_stream;
      string_stream << std::scientific << std::setprecision(precision)
                    << to_convert.second.at(i).real();
      vector_entry.replace(index, real_key_string.size(), string_stream.str());
    }
    if (std::size_t index = vector_entry.find(imag_key_string);
        index != std::string::npos) {
      std::ostringstream string_stream;
      string_stream << std::scientific << std::setprecision(precision)
                    << to_convert.second.at(i).imag();
      vector_entry.replace(index, imag_key_string.size(), string_stream.str());
    }
    full_vector_string += vector_entry;
    if (i + 1 != static_cast<int>(to_convert.second.size()))
      full_vector_string += delimiter;
  }

  std::string vector_string = output_term_to_string_map.at(OutputTerm::kVector);
  if (std::size_t index = output.find(vector_string);
      index != std::string::npos) {
    output.replace(index, vector_string.size(), full_vector_string);
  }

  return output;
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_converter_ptr;
  EXPECT_NO_THROW({
    test_converter_ptr = std::make_unique<ConverterType>();
  });
  EXPECT_EQ(test_converter_ptr->output_format(), default_output_);
  EXPECT_EQ(test_converter_ptr->precision(), default_precision_);
  EXPECT_EQ(test_converter_ptr->vector_delimiter(), default_delimiter_);
  EXPECT_EQ(test_converter_ptr->vector_output_format(),
            default_vector_entry_output_format);
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest, SetPrecision) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.precision(), default_precision_);
  const int new_precision(test_helpers::RandomInt(3,5));
  test_converter.set_precision(new_precision);
  EXPECT_EQ(test_converter.precision(), new_precision);
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest, SetDelimiter) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.vector_delimiter(), default_delimiter_);
  const std::string new_delimiter("; ");
  test_converter.set_vector_delimiter(new_delimiter);
  EXPECT_EQ(test_converter.vector_delimiter(), new_delimiter);
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
       SetVectorOutputFormat) {
  ConverterType test_converter;
  std::string new_format{"Entry ${INDEX}, Real: ${REAL}, ${IMAG}"};
  auto returned_string = test_converter.SetVectorOutputFormat({
        "Entry ", VectorTerm::kIndex, ", Real: ", VectorTerm::kReal,
        ", ", VectorTerm::kImag});
  EXPECT_EQ(returned_string, new_format);
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
    ConvertDefaults) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.Convert(test_input),
            GetExpectedOutput(test_input,
                              test_converter,
                              default_output_,
                              default_vector_entry_output_format,
                              default_delimiter_,
                              default_precision_));
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
       ConvertDifferentOutputFormat) {
  ConverterType test_converter;
  std::string new_format{"Vector: ${VECTOR} is paired with int ${INT}\n"};
  test_converter.SetOutputFormat(
      {"Vector: ", OutputTerm::kVector, " is paired with int ", OutputTerm::kInt,
       "\n"});

  EXPECT_EQ(test_converter.Convert(test_input),
            GetExpectedOutput(test_input,
                              test_converter,
                              new_format,
                              default_vector_entry_output_format,
                              default_delimiter_,
                              default_precision_));
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
       ConvertDifferentVectorOutputFormat) {
  ConverterType test_converter;
  std::string new_format{"[Entry ${INDEX}: Real: ${REAL}, Imag: ${IMAG}]"};
  test_converter.SetVectorOutputFormat({
    "[Entry ", VectorTerm::kIndex, ": Real: ", VectorTerm::kReal,
    ", Imag: ", VectorTerm::kImag, "]"});

  EXPECT_EQ(test_converter.Convert(test_input),
            GetExpectedOutput(test_input,
                              test_converter,
                              default_output_,
                              new_format,
                              default_delimiter_,
                              default_precision_));
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
       ConvertDifferentDelimiter) {
  ConverterType test_converter;
  std::string new_delimiter{";"};
  test_converter.set_vector_delimiter(new_delimiter);
  EXPECT_EQ(test_converter.Convert(test_input),
            GetExpectedOutput(test_input,
                              test_converter,
                              default_output_,
                              default_vector_entry_output_format,
                              new_delimiter,
                              default_precision_));
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest,
       ConverterDifferentPrecision) {
  ConverterType test_converter;
  const int new_precision{test_helpers::RandomInt(5, 10)};
  test_converter.set_precision(new_precision);
  EXPECT_EQ(test_converter.Convert(test_input),
            GetExpectedOutput(test_input,
                              test_converter,
                              default_output_,
                              default_vector_entry_output_format,
                              default_delimiter_,
                              new_precision));
}




} // namespace
