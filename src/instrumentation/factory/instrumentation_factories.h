#ifndef BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_

#include <memory>

#include "instrumentation/instrument.h"
#include "instrumentation/converter/converter_i.h"
#include "instrumentation/outstream/outstream_i.h"

namespace bart {

namespace instrumentation {

namespace factory {
template <typename InputType, typename OutputType>
using ConverterType = instrumentation::converter::ConverterI<InputType, OutputType>;
template <typename InputType>
using InstrumentType = instrumentation::InstrumentI<InputType>;
template <typename OutputType>
using OutstreamType = instrumentation::outstream::OutstreamI<OutputType>;

template <typename InputType>
std::unique_ptr<InstrumentType<InputType>> MakeBasicInstrument(
    std::unique_ptr<OutstreamType<InputType>> outstream_ptr_);

template <typename InputType, typename OutputType, typename ... DependencyTypes>
std::unique_ptr<ConverterType<InputType, OutputType>> MakeConverter(
    DependencyTypes ... dependencies);

template <typename InputType, typename OutputType>
std::unique_ptr<InstrumentType<InputType>> MakeInstrument(
    std::unique_ptr<ConverterType<InputType, OutputType>> converter_ptr_,
    std::unique_ptr<OutstreamType<OutputType>> outstream_ptr_);

template <typename OutputType, typename ... DependencyTypes>
std::unique_ptr<OutstreamType<OutputType>> MakeOutstream(
    DependencyTypes ... dependencies);


} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
