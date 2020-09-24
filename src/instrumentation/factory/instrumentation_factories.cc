#include "instrumentation_factories.h"

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>
#include <deal.II/base/conditional_ostream.h>

#include "calculator/fourier/fourier_transform_i.h"
#include "convergence/status.h"
#include "instrumentation/basic_instrument.h"
#include "instrumentation/instrument.h"
#include "instrumentation/converter/dealii_to_complex_vector.h"
#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/to_string/convergence_to_string.h"
#include "instrumentation/converter/to_string/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/to_string/int_double_pair_to_string.h"
#include "instrumentation/converter/to_string/string_color_pair_to_string.h"
#include "instrumentation/outstream/to_conditional_ostream.h"
#include "instrumentation/outstream/to_ostream.h"
#include "utility/colors.h"

namespace bart {

namespace instrumentation {

namespace factory {


template<typename InputType>
std::unique_ptr<InstrumentType<InputType>> MakeBasicInstrument(std::unique_ptr<
    OutstreamType<InputType>> outstream_ptr_) {
  using ReturnType = instrumentation::BasicInstrument<InputType>;
  return std::make_unique<ReturnType>(std::move(outstream_ptr_));
}

// MAKE CONVERTER ==============================================================

namespace  {
  // types used by converters
using ComplexVector = std::vector<std::complex<double>>;
using DealiiVector = dealii::Vector<double>;
using IntComplexVectorPair = std::pair<int, std::vector<std::complex<double>>>;
using IntDoublePair = std::pair<int, double>;
using StringColorPair = std::pair<std::string, utility::Color>;

template <typename InputType>
using ConvertThisToStringPtr = std::unique_ptr<ConverterType<InputType, std::string>>;

// Dependencies
using FourierCalculatorPtr = std::unique_ptr<calculator::fourier::FourierTransformI>;
} // namespace

template <>
std::unique_ptr<ConverterType<DealiiVector, ComplexVector>>
MakeConverter<DealiiVector, ComplexVector>() {
  using ReturnType = instrumentation::converter::DealiiToComplexVector;
  return std::make_unique<ReturnType>();
}

// FourierConverters ===========================================================

template <>
std::unique_ptr<ConverterType<ComplexVector, ComplexVector>>
MakeConverter<ComplexVector, ComplexVector, FourierCalculatorPtr>(
    FourierCalculatorPtr fourier_calculator_ptr) {
using ReturnType = instrumentation::converter::fourier::FourierTransform;
  return std::make_unique<ReturnType>(std::move(fourier_calculator_ptr));
}

// ToStringConverters ==========================================================

template <>
std::unique_ptr<ConverterType<convergence::Status, std::string>>
MakeConverter<convergence::Status, std::string>() {
  using ReturnType = instrumentation::converter::to_string::ConvergenceToString;
  return std::make_unique<ReturnType>();
}

template <>
auto MakeConverter<IntComplexVectorPair, std::string, int> (const int precision)
-> ConvertThisToStringPtr<IntComplexVectorPair> {
    using ReturnType = instrumentation::converter::to_string::IntVectorComplexPairToString;
    auto return_ptr = std::make_unique<ReturnType>();
    return_ptr->set_precision(precision);
    return return_ptr;
}

template <>
ConvertThisToStringPtr<IntDoublePair>
MakeConverter<IntDoublePair, std::string, int>(const int precision) {
  using ReturnType = instrumentation::converter::to_string::IntDoublePairToString;
  auto return_ptr = std::make_unique<ReturnType>();
  return_ptr->set_precision(precision);
  return return_ptr;
}

template <>
ConvertThisToStringPtr<StringColorPair>
MakeConverter<StringColorPair, std::string>() {
  using ReturnType = instrumentation::converter::to_string::StringColorPairToString;
  auto return_ptr = std::make_unique<ReturnType>();
  return return_ptr;
}

// MAKE INSTRUMENT =============================================================

template <typename InputType, typename OutputType>
std::unique_ptr<InstrumentType<InputType>> MakeInstrument(
    std::unique_ptr<ConverterType<InputType, OutputType>> converter_ptr_,
    std::unique_ptr<OutstreamType<OutputType>> outstream_ptr_) {
  using ReturnType = instrumentation::Instrument<InputType, OutputType>;
  return std::make_unique<ReturnType>(std::move(converter_ptr_),
                                      std::move(outstream_ptr_));
}

// MAKE OUTSTREAM ==============================================================

template <>
std::unique_ptr<OutstreamType<std::string>>
MakeOutstream<std::string, std::unique_ptr<dealii::ConditionalOStream>>(
    std::unique_ptr<dealii::ConditionalOStream> conditional_ostream_ptr) {
  using ReturnType = instrumentation::outstream::ToConditionalOstream;
  return std::make_unique<ReturnType>(std::move(conditional_ostream_ptr));
}

template <>
std::unique_ptr<OutstreamType<std::string>>
MakeOutstream<std::string, std::unique_ptr<std::ostream>>(
    std::unique_ptr<std::ostream> ostream_ptr) {
  using ReturnType = instrumentation::outstream::ToOstream;
  return std::make_unique<ReturnType>(std::move(ostream_ptr));
}

template std::unique_ptr<InstrumentType<convergence::Status>> MakeInstrument(std::unique_ptr<ConverterType<convergence::Status, std::string>>, std::unique_ptr<OutstreamType<std::string>>);
template std::unique_ptr<InstrumentType<std::pair<int, double>>> MakeInstrument(std::unique_ptr<ConverterType<std::pair<int, double>, std::string>>, std::unique_ptr<OutstreamType<std::string>>);
template std::unique_ptr<InstrumentType<std::pair<std::string, utility::Color>>> MakeInstrument(std::unique_ptr<ConverterType<std::pair<std::string, utility::Color>, std::string>>, std::unique_ptr<OutstreamType<std::string>>);
template std::unique_ptr<InstrumentType<std::string>> MakeBasicInstrument(std::unique_ptr<OutstreamType<std::string>>);

} // namespace factory

} // namespace instrumentation

} // namespace bart
