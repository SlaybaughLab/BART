#include "instrumentation/factory/component_factories.h"

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>
#include <deal.II/base/conditional_ostream.h>

#include "calculator/fourier/tests/fourier_transform_mock.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/to_string/convergence_to_string.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/to_string/int_double_pair_to_string.h"
#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/converter/tests/converter_mock.h"
#include "instrumentation/instrument.h"
#include "instrumentation/basic_instrument.h"
#include "instrumentation/outstream/to_ostream.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/outstream/tests/outstream_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/colors.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace convergence = bart::convergence;
namespace utility = bart::utility;
using ::testing::WhenDynamicCastTo, ::testing::NotNull;

class ComponentFactoriesIntegrationTests : public ::testing::Test {
 public:
};

TEST_F(ComponentFactoriesIntegrationTests, ConditionalOstream) {
  auto outstream_ptr = instrumentation::factory::MakeOutstream<std::string>(
      std::make_unique<dealii::ConditionalOStream>(std::cout, false));
  ASSERT_NE(outstream_ptr, nullptr);
  using ExpectedType = instrumentation::outstream::ToConditionalOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

TEST_F(ComponentFactoriesIntegrationTests, ToOStreamOutStream) {
  using namespace instrumentation::factory;
  auto outstream_ptr = MakeOutstream<std::string, std::unique_ptr<std::ostream>>(
      std::make_unique<std::ostringstream>());
  ASSERT_NE(outstream_ptr, nullptr);
  using ExpectedType = instrumentation::outstream::ToOstream;
  ASSERT_NE(dynamic_cast<ExpectedType*>(outstream_ptr.get()), nullptr);
}

// CONVERTER TESTS =============================================================

TEST_F(ComponentFactoriesIntegrationTests, ConverterDealiiVectorToComplexVector) {
  using ComplexVector = std::vector<std::complex<double>>;
  using DealiiVector = dealii::Vector<double>;
  using ExpectedType = instrumentation::converter::DealiiToComplexVector;

  auto converter_ptr = instrumentation::factory::MakeConverter<DealiiVector, ComplexVector>();
  ASSERT_NE(converter_ptr, nullptr);
  EXPECT_THAT(converter_ptr.get(), WhenDynamicCastTo<ExpectedType*>(NotNull()));
}

// FOURIER CONVERTERS ==========================================================

TEST_F(ComponentFactoriesIntegrationTests, ConverterFourier) {
  using ComplexVector = std::vector<std::complex<double>>;
  using FourierCalculator = bart::calculator::fourier::FourierTransformMock;
  using FourierCalculatorI = bart::calculator::fourier::FourierTransformI;
  using ExpectedType = instrumentation::converter::fourier::FourierTransform;
  std::unique_ptr<FourierCalculatorI> fourier_calculator_ptr =
      std::make_unique<FourierCalculator>();

  auto converter_ptr =
      instrumentation::factory::MakeConverter<ComplexVector, ComplexVector>(
          std::move(fourier_calculator_ptr));
  ASSERT_NE(converter_ptr, nullptr);
  EXPECT_THAT(converter_ptr.get(), WhenDynamicCastTo<ExpectedType*>(NotNull()));
}

// TO_STRING CONVERTERS ========================================================

TEST_F(ComponentFactoriesIntegrationTests, ConverterConvergenceToString) {
  auto converter_ptr =  instrumentation::factory::MakeConverter<convergence::Status, std::string>();
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::to_string::ConvergenceToString;
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}

TEST_F(ComponentFactoriesIntegrationTests,
    ConverterIntComplexVectorPairToString) {
  const int precision = 5;
  using InputType = std::pair<int, std::vector<std::complex<double>>>;
  using OutputType = std::string;
  using ExpectedReturnType = instrumentation::converter::to_string::IntVectorComplexPairToString;
  auto converter_ptr = instrumentation::factory::MakeConverter<InputType, OutputType>(precision);
  auto dynamic_ptr = dynamic_cast<ExpectedReturnType*>(converter_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  EXPECT_EQ(dynamic_ptr->precision(), 5);
}

TEST_F(ComponentFactoriesIntegrationTests, ConverterIntDoubleToString) {
  const int precision = 5;
  auto converter_ptr = instrumentation::factory::MakeConverter<std::pair<int, double>, std::string>(precision);
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::to_string::IntDoublePairToString;
  auto dynamic_ptr = dynamic_cast<ExpectedType*>(converter_ptr.get());
  ASSERT_NE(dynamic_ptr, nullptr);
  EXPECT_EQ(dynamic_ptr->precision(), precision);
}

TEST_F(ComponentFactoriesIntegrationTests, StringColorPairToString) {
  using StringColorPair = std::pair<std::string, utility::Color>;
  auto converter_ptr =  instrumentation::factory::MakeConverter<StringColorPair,
                                                                std::string>();
  ASSERT_NE(converter_ptr, nullptr);
  using ExpectedType = instrumentation::converter::to_string::StringColorPairToString;
  ASSERT_NE(dynamic_cast<ExpectedType*>(converter_ptr.get()), nullptr);
}


TEST_F(ComponentFactoriesIntegrationTests, MakeBasicInstrument) {
  using InputType = std::string;
  using OutstreamType = instrumentation::outstream::OutstreamMock<InputType>;

  auto instrument_ptr = instrumentation::factory::MakeBasicInstrument<InputType>(
      std::make_unique<OutstreamType>());
  using ExpectedType = instrumentation::BasicInstrument<InputType>;
  ASSERT_NE(instrument_ptr, nullptr);
  ASSERT_NE(dynamic_cast<ExpectedType*>(instrument_ptr.get()), nullptr);
}

TEST_F(ComponentFactoriesIntegrationTests, MakeInstrument) {
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
