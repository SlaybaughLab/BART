#ifndef BART_SRC_INSTRUMENTATION_PORT_H_
#define BART_SRC_INSTRUMENTATION_PORT_H_

#include "instrumentation/instrument_i.h"

#include <memory>

namespace bart {

namespace instrumentation {

template <typename DataType, typename DataName>
class Port {
 public:
  using InstrumentType = InstrumentI<DataType>;
  virtual ~Port() = default;
  void AddInstrument(std::shared_ptr<InstrumentType> instrument_ptr) {
    instrument_ptr_ = instrument_ptr; }

  InstrumentType* instrument_ptr() { return instrument_ptr_.get(); }
 protected:
  std::shared_ptr<InstrumentType> instrument_ptr_ = nullptr;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_PORT_H_
