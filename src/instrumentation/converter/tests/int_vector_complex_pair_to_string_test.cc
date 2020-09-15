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

  const int default_precision_{2};
  const std::string default_output_{"${INT}, ${VECTOR}\n"};
  const std::string default_vector_entry_output_format{"(${REAL}, ${IMAG}$)"};
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

} // namespace
