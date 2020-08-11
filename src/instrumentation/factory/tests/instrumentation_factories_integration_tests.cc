#include "instrumentation/factory/instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "instrumentation/converter/convergence_to_string.h"
#include "instrumentation/converter/int_double_pair_to_string.h"
#include "instrumentation/converter/tests/converter_mock.h"
#include "instrumentation/instrument.h"
#include "instrumentation/basic_instrument.h"
#include "instrumentation/outstream/to_ostream.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/outstream/tests/outstream_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace convergence = bart::convergence;

class InstrumentationFactoriesIntegrationTests : public ::testing::Test {
 public:
};

TEST_F(InstrumentationFactoriesIntegrationTests, ConditionalOstream) {
  auto outstream_ptr = instrumentation::factory::MakeOutstream<std::string>(
      std::make_unique<dealii::ConditionalOStream>(std::cout, false));
  ASSERT_NE(outstream_ptr, nullptr);
  using ExpectedType = instrumentation::outstream::ToConditionalOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

TEST_F(InstrumentationFactoriesIntegrationTests, ToOStreamOutStream) {
  using namespace instrumentation::factory;
  auto outstream_ptr = MakeOutstream<std::string, std::unique_ptr<std::ostream>>(
      std::make_unique<std::ostringstream>());
  ASSERT_NE(outstream_ptr, nullptr);
  using ExpectedType = instrumentation::outstream::ToOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}



TEST_F(InstrumentationFactoriesIntegrationTests, ConverterConvergenceToString) {
  auto converter_ptr =  instrumentation::factory::MakeConverter<convergence::Status, std::string>();
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::ConvergenceToString;
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(InstrumentationFactoriesIntegrationTests, ConverterIntDoubleToString) {
  const int precision = 5;
  auto converter_ptr = instrumentation::factory::MakeConverter<std::pair<int, double>, std::string>(precision);
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::IntDoublePairToString;
  auto dynamic_ptr = dynamic_cast<ExpectedType*>(converter_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  EXPECT_EQ(dynamic_ptr->precision(), precision);
}

TEST_F(InstrumentationFactoriesIntegrationTests, MakeBasicInstrument) {
  using InputType = std::string;
  using OutstreamType = instrumentation::outstream::OutstreamMock<InputType>;

  auto instrument_ptr = instrumentation::factory::MakeBasicInstrument<InputType>(
      std::make_unique<OutstreamType>());
  using ExpectedType = instrumentation::BasicInstrument<InputType>;
  ASSERT_NE(instrument_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(instrument_ptr.get()), nullptr);
}

TEST_F(InstrumentationFactoriesIntegrationTests, MakeInstrument) {
  using InputType = convergence::Status;
  using OutputType = std::string;
  using ConverterType = instrumentation::converter::ConverterMock<InputType, OutputType>;
  using OutstreamType = instrumentation::outstream::OutstreamMock<OutputType>;

  auto instrument_ptr =
      instrumentation::factory::MakeInstrument<InputType, OutputType>(
          std::make_unique<ConverterType>(), std::make_unique<OutstreamType>());
  using ExpectedType = instrumentation::Instrument<InputType, OutputType>;
  ASSERT_NE(instrument_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(instrument_ptr.get()), nullptr);
}

} // namespace
