#include "instrumentation/converter/string_color_pair_to_string.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

class InstrumentationConverterStringColorPairToStringTest
 : public ::testing::Test {
 public:
  using ConverterType = instrumentation::converter::StringColorPairToString;

  std::string default_output{"${COLOR_CODE}${STRING}${COLOR_RESET_CODE}"};
};

TEST_F(InstrumentationConverterStringColorPairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_ptr_;
  EXPECT_NO_THROW({test_ptr_ = std::make_unique<ConverterType>();});
  EXPECT_EQ(test_ptr_->output_format(), default_output);
}

} // namespace
