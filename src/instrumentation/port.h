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
  void Expose(const DataType& data_to_expose) {
    if (instrument_ptr_ != nullptr) {
      instrument_ptr_->Read(data_to_expose);
    }
  }
  InstrumentType* instrument_ptr() { return instrument_ptr_.get(); }
  struct Use_dot_after_GetPort_not_arrow {};
 protected:
  std::shared_ptr<InstrumentType> instrument_ptr_ = nullptr;
 private:
  Use_dot_after_GetPort_not_arrow* operator->() { return nullptr; }
};

template <typename PortType, typename T>
PortType& GetPort(T& to_expose) {
  static_assert(std::is_polymorphic_v<T>,
      "Error in GetPort, passed object must be a reference to a polymorphic "
      "class that can be cast to a Port. Pointers cannot be passed, you "
      "must dereference to the underlying object using *");
  return dynamic_cast<PortType&>(to_expose);
}

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_PORT_H_
