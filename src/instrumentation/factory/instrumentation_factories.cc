#include "instrumentation_factories.h"

#include <deal.II/base/conditional_ostream.h>

#include "convergence/status.h"
#include "instrumentation/basic_instrument.h"
#include "instrumentation/instrument.h"
#include "instrumentation/converter/convergence_to_string.h"
#include "instrumentation/outstream/to_conditional_ostream.h"

namespace bart {

namespace instrumentation {

namespace factory {

template<typename InputType>
std::unique_ptr<InstrumentType<InputType>> MakeBasicInstrument(std::unique_ptr<
    OutstreamType<InputType>> outstream_ptr_) {
  using ReturnType = instrumentation::BasicInstrument<InputType>;
  return std::make_unique<ReturnType>(std::move(outstream_ptr_));
}

template <>
std::unique_ptr<ConverterType<convergence::Status, std::string>>
MakeConverter<convergence::Status, std::string>() {
  using ReturnType = instrumentation::converter::ConvergenceToString;
  return std::make_unique<ReturnType>();
}

template <typename InputType, typename OutputType>
std::unique_ptr<InstrumentType<InputType>> MakeInstrument(
    std::unique_ptr<ConverterType<InputType, OutputType>> converter_ptr_,
    std::unique_ptr<OutstreamType<OutputType>> outstream_ptr_) {
  using ReturnType = instrumentation::Instrument<InputType, OutputType>;
  return std::make_unique<ReturnType>(std::move(converter_ptr_),
                                      std::move(outstream_ptr_));
}

template <>
std::unique_ptr<OutstreamType<std::string>>
MakeOutstream<std::string, std::unique_ptr<dealii::ConditionalOStream>>(
    std::unique_ptr<dealii::ConditionalOStream> conditional_ostream_ptr) {
  using ReturnType = instrumentation::outstream::ToConditionalOstream;
  return std::make_unique<ReturnType>(std::move(conditional_ostream_ptr));
}

template std::unique_ptr<InstrumentType<convergence::Status>> MakeInstrument(std::unique_ptr<ConverterType<convergence::Status, std::string>>, std::unique_ptr<OutstreamType<std::string>>);
template std::unique_ptr<InstrumentType<std::string>> MakeBasicInstrument(std::unique_ptr<OutstreamType<std::string>>);

} // namespace factory

} // namespace instrumentation

} // namespace bart
