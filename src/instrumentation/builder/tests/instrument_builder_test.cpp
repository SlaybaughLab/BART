#include "instrumentation/builder/instrument_builder.hpp"

#include <string>
#include <utility>

#include "convergence/status.h"
#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/converter/to_string/convergence_to_string.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/instrument.h"
#include "instrumentation/basic_instrument.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

namespace instrumentation = bart::instrumentation;

class InstrumentationBuilderInstrumentBuilderTest : public ::testing::Test {
 public:
  using Builder = instrumentation::builder::InstrumentBuilder;
  using InstrumentName = instrumentation::builder::InstrumentName;
  using ConvergenceStatus = bart::convergence::Status;
  using StringColorPair = std::pair<std::string, bart::utility::Color>;
};

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       ColorStringToConditionalOstream) {
  using InstrumentType = instrumentation::Instrument<StringColorPair, std::string>;
  using ConverterType = instrumentation::converter::to_string::StringColorPairToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       ConvergenceStatusToConditionalOstream) {
  using InstrumentType = instrumentation::Instrument<ConvergenceStatus, std::string>;
  using ConverterType = instrumentation::converter::to_string::ConvergenceToString;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<ConvergenceStatus>(
      InstrumentName::kConvergenceStatusToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ConverterType*>(dynamic_ptr->converter_ptr()), nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       StringToConditionalOstream) {
  using InstrumentType = instrumentation::BasicInstrument<std::string>;
  using OutstreamType = instrumentation::outstream::ToConditionalOstream;
  auto instrument_ptr = Builder::BuildInstrument<std::string>(
      instrumentation::builder::InstrumentName::kStringToConditionalOstream);
  ASSERT_NE(instrument_ptr, nullptr);
  auto dynamic_ptr = dynamic_cast<InstrumentType*>(instrument_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  ASSERT_NE(dynamic_cast<OutstreamType*>(dynamic_ptr->outstream_ptr()), nullptr);
}

TEST_F(InstrumentationBuilderInstrumentBuilderTest,
       BadInstrumentNamesPerType) {
  EXPECT_ANY_THROW(Builder::BuildInstrument<std::string>(
      instrumentation::builder::InstrumentName::kColorStatusToConditionalOstream));
  EXPECT_ANY_THROW(Builder::BuildInstrument<StringColorPair>(
      instrumentation::builder::InstrumentName::kStringToConditionalOstream));

  for (auto bad_name : {InstrumentName::kColorStatusToConditionalOstream,
                        InstrumentName::kStringToConditionalOstream}) {
    EXPECT_ANY_THROW(Builder::BuildInstrument<ConvergenceStatus>(bad_name));
  }

}



} // namespace
