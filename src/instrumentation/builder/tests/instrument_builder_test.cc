#include "instrumentation/builder/instrument_builder.h"

#include <string>
#include <utility>

#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/instrument.h"

#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationBuilderInstrumentBuilderTest : public ::testing::Test {
 public:
  instrumentation::builder::InstrumentBuilder builder_;
};

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       ColorStringToConditionalOstream) {
  using StringColorPair = std::pair<std::string, bart::utility::Color>;
  using InstrumentType = instrumentation::Instrument<StringColorPair, std::string>;
  using ConverterType = instrumentation::converter::to_string::StringColorPairToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = builder_.BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);

}

} // namespace
