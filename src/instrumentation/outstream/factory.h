#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_H_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_H_

#include "utility/factory/auto_registering_factory.h"

namespace bart::instrumentation::outstream {

template <typename DataType> class OutstreamI;

enum class OutstreamName {
  kToOstream,
  kToConditionalOstream,
};

template <typename DataType, typename...ArgTypes>
class OutstreamIFactory : public utility::factory::AutoRegisteringFactory<
    OutstreamName,
    std::unique_ptr<OutstreamI<DataType>>(*)(ArgTypes...)> {};

} // namespace bart::instrumentation::outstream

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_H_
