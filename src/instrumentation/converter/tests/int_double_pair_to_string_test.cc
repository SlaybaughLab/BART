#include "instrumentation/converter/int_double_pair_to_string.h"

#include <optional>

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
//  std::string GetExpectedOutput(const std::pair<int, double>, ConverterType&,
//                                const std::string, const int) const;
};

//std::string
//InstrumentationConverterIntDoublePairToStringTest::GetExpectedOutput(
//    const std::pair<int, double> /*to_convert*/,
//    ConverterType& /*converter*/,
//    const std::string /*format*/,
//    const int /*precision*/) const {
//  using OutputTerm = ConverterType::OutputTerm;
//}

TEST_F(InstrumentationConverterIntDoublePairToStringTest, Constructor) {
  std::unique_ptr<ConverterType> test_converter_ptr;
  EXPECT_NO_THROW({
    test_converter_ptr = std::make_unique<ConverterType>();
  });
  EXPECT_EQ(test_converter_ptr->output_format(), "${INDEX}, ${VALUE}\n");
  EXPECT_EQ(test_converter_ptr->precision(), 2);
}

} // namespace
