#ifndef BART_SRC_INSTRUMENTATION_INSTRUMENT_ARRAY_HPP_
#define BART_SRC_INSTRUMENTATION_INSTRUMENT_ARRAY_HPP_

#include "instrumentation/instrument_i.h"
#include "utility/uncopyable.h"

namespace bart::instrumentation {


/*! \brief An instrument that can hold multiple other instruments to call.
 *
 * Derives from the base instrument, so it can be used anywhere that a normal
 * instrument may be called. This is best used for reading out to multiple
 * outstreams or performing different conversions on the same data.
 *
 * @tparam InputType type of data to be read.
 */
template <typename InputType>
class InstrumentArray : public InstrumentI<InputType>,
                        public utility::Uncopyable {
 public:
  void Read(const InputType &input) override;
};


} // namespace bart::instrumentation

#endif //BART_SRC_INSTRUMENTATION_INSTRUMENT_ARRAY_HPP_
