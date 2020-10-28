#ifndef BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_HPP_
#define BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_HPP_

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

[[nodiscard]] inline auto to_string(OutstreamName to_convert) -> std::string {
  switch (to_convert) {
    case (OutstreamName::kToOstream):
      return std::string{"OutstreamName::kToOstream"};
    case (OutstreamName::kToConditionalOstream):
      return std::string{"OutstreamName::kToConditionalOstream"};
  }
  return std::string{"Unknown OutstreamName conversion to string requested."};
}

} // namespace bart::instrumentation::outstream

#endif //BART_SRC_INSTRUMENTATION_OUTSTREAM_FACTORY_HPP_
