#ifndef BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_
#define BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_

#include <memory>

#include "instrumentation/instrument_i.h"

namespace bart::instrumentation::builder {

enum class InstrumentName {
  kColorStatusToConditionalOstream = 0,
  kConvergenceStatusToConditionalOstream = 1,
  kFourierTransformOfSingleGroupScalarFluxErrorToFile = 2,
  kFourierTransformOfAllGroupScalarFluxErrorToFile = 3,
  kIntDoublePairToFile = 4,
  kStringToConditionalOstream = 5,
  kDoubleToFile = 6,
};

class InstrumentBuilder {
 public:
  template <typename InputType, typename ...Args>
  [[nodiscard]] static auto BuildInstrument(
      const InstrumentName name, Args ... args) -> std::unique_ptr<InstrumentI<InputType>>;
};

} // namespace bart::instrumentation::builder

#endif //BART_SRC_INSTRUMENTATION_BUILDER_INSTRUMENT_BUILDER_HPP_
