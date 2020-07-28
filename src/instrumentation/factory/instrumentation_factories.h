#ifndef BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
#define BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_

#include <memory>

#include "instrumentation/output/output_i.h"

namespace bart {

namespace instrumentation {

namespace factory {

template <typename OutputType, typename ... DependencyTypes>
std::unique_ptr<instrumentation::output::OutputI<OutputType>> MakeOutputter(
    DependencyTypes ... dependencies);

} // namespace factory

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_FACTORY_INSTRUMENTATION_FACTORIES_H_
