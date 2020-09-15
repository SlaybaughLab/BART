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
}

} // namespace
