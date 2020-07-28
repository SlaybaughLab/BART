#ifndef BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_

#include <memory>

#include "instrumentation/converter/converter_i.h"
#include "instrumentation/outstream/outstream_i.h"

namespace bart {

namespace instrumentation {

namespace factory {

template <typename OutputType>
using OutstreamType = instrumentation::outstream::OutstreamI<OutputType>;

template <typename InputType, typename OutputType>
using ConverterType = instrumentation::converter::ConverterI<InputType, OutputType>;

template <typename OutputType, typename ... DependencyTypes>
std::unique_ptr<OutstreamType<OutputType>> MakeOutstream(
    DependencyTypes ... dependencies);

template <typename InputType, typename OutputType, typename ... DependencyTypes>
std::unique_ptr<ConverterType<InputType, OutputType>> MakeConverter(
    DependencyTypes ... dependencies);

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
