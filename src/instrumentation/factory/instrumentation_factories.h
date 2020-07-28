#ifndef BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_

#include <memory>

#include "instrumentation/outstream/outstream_i.h"

namespace bart {

namespace instrumentation {

namespace factory {

template <typename OutputType>
using OutStreamType = instrumentation::outstream::OutstreamI<OutputType>;

template <typename OutputType, typename ... DependencyTypes>
std::unique_ptr<OutStreamType<OutputType>> MakeOutstream(
    DependencyTypes ... dependencies);

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
