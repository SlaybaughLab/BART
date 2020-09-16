#include "instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "convergence/status.h"
#include "instrumentation/basic_instrument.h"
#include "instrumentation/instrument.h"
#include "instrumentation/converter/convergence_to_string.h"
#include "instrumentation/converter/int_vector_complex_pair_to_string.h"
#include "instrumentation/converter/int_double_pair_to_string.h"
#include "instrumentation/converter/string_color_pair_to_string.h"
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

template <>
std::unique_ptr<ConverterType<convergence::Status, std::string>>
MakeConverter<convergence::Status, std::string>() {
  using ReturnType = instrumentation::converter::ConvergenceToString;
  return std::make_unique<ReturnType>();
}

template <>
std::unique_ptr<ConverterType<std::pair<int, double>, std::string>>
MakeConverter<std::pair<int, double>, std::string, int>(const int precision) {
  using ReturnType = instrumentation::converter::IntDoublePairToString;
  auto return_ptr = std::make_unique<ReturnType>();
  return_ptr->set_precision(precision);
  return return_ptr;
}

template <>
std::unique_ptr<ConverterType<std::pair<std::string, utility::Color>, std::string>>
MakeConverter<std::pair<std::string, utility::Color>, std::string>() {
  using ReturnType = instrumentation::converter::StringColorPairToString;
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
