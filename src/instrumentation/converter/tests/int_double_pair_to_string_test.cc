#include "instrumentation/converter/int_double_pair_to_string.h"

#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterIntDoublePairToStringTest
    : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::IntDoublePairToString;
};

TEST_F(InstrumentationConverterIntDoublePairToStringTest, Constructor) {
  EXPECT_NO_THROW({
    ConverterType new_converter;
  });
}

} // namespace
