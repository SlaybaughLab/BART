#include "instrumentation/instrument.h"

#include <string>

#include "convergence/status.h"

namespace bart {

namespace instrumentation {

template<typename InputType, typename OutputType>
Instrument<InputType, OutputType>::Instrument(
    std::unique_ptr<ConverterType> converter_ptr,
    std::unique_ptr<OutputterType> outputter_ptr)
    : converter_ptr_(std::move(converter_ptr)),
      outputter_ptr_(std::move(outputter_ptr)) {}

template class Instrument<convergence::Status, std::string>;

} // namespace instrumentation

} // namespace bart
