#ifndef BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_
#define BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_

#include <memory>

#include "instrumentation/instrument_i.h"

namespace bart::instrumentation::builder {

enum class InstrumentName {
  kColorStatusToConditionalOstream = 0,
  kStringToConditionalOstream = 1,
};

class InstrumentBuilder {
 public:
  template <typename InputType>
  [[nodiscard]] std::unique_ptr<InstrumentI<InputType>> BuildInstrument(
      InstrumentName name) const;
};

} // namespace bart::instrumentation::builder

#endif //BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_
