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

  const int default_precision_{2};
  const std::string default_output_{"${INT}, (${REAL}, ${IMAG})"};
};

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_converter_ptr;
  EXPECT_NO_THROW({
    test_converter_ptr = std::make_unique<ConverterType>();
  });
  EXPECT_EQ(test_converter_ptr->output_format(), default_output_);
  EXPECT_EQ(test_converter_ptr->precision(), default_precision_);
}

TEST_F(InstrumentationConverterIntVectorComplexPairToStringTest, SetPrecision) {
  ConverterType test_converter;
  EXPECT_EQ(test_converter.precision(), default_precision_);
  const int new_precision(test_helpers::RandomInt(3,5));
  test_converter.set_precision(new_precision);
  EXPECT_EQ(test_converter.precision(), new_precision);
}

} // namespace
