#ifndef BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
#define BART_SRC_INSTRUMENTATION_INSTRUMENT_H_

#include "instrumentation/instrument_i.h"

namespace bart {

namespace instrumentation {

template <typename InputType, typename OutputType>
class Instrument : public InstrumentI<InputType> {
 public:
  void Read(const InputType &input) override {}
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
