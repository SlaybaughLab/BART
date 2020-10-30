#ifndef BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
#define BART_SRC_INSTRUMENTATION_INSTRUMENT_H_

#include <memory>

#include "instrumentation/instrument_i.h"
#include "instrumentation/converter/converter_i.h"
#include "instrumentation/outstream/outstream_i.h"


namespace bart {

namespace instrumentation {

template <typename InputType, typename OutputType>
class Instrument : public InstrumentI<InputType> {
 public:
  using ConverterType = converter::ConverterI<InputType, OutputType>;
  using OutstreamType = outstream::OutstreamI<OutputType>;

  Instrument(std::unique_ptr<ConverterType>, std::unique_ptr<OutstreamType>);
  virtual ~Instrument() = default;

  virtual void Read(const InputType &input) override;

  ConverterType* converter_ptr() { return converter_ptr_.get(); }
  OutstreamType* outstream_ptr() { return outstream_ptr_.get(); }
 protected:
  template <typename T>
  void AssertNotNull(T* ptr, std::string, std::string);
  std::unique_ptr<ConverterType> converter_ptr_ = nullptr;
  std::unique_ptr<OutstreamType> outstream_ptr_ = nullptr;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
