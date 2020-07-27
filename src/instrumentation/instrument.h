#ifndef BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
#define BART_SRC_INSTRUMENTATION_INSTRUMENT_H_

#include <memory>

#include "instrumentation/instrument_i.h"
#include "instrumentation/converter/converter_i.h"
#include "instrumentation/output/output_i.h"


namespace bart {

namespace instrumentation {

template <typename InputType, typename OutputType>
class Instrument : public InstrumentI<InputType> {
 public:
  using ConverterType = converter::ConverterI<InputType, OutputType>;
  using OutputterType = output::OutputI<OutputType>;

  Instrument(std::unique_ptr<ConverterType>, std::unique_ptr<OutputterType>);
  void Read(const InputType &input) override {}

  ConverterType* converter_ptr() { return converter_ptr_.get(); }
  OutputterType* outputter_ptr() { return outputter_ptr_.get(); }
 private:
  std::unique_ptr<ConverterType> converter_ptr_ = nullptr;
  std::unique_ptr<OutputterType> outputter_ptr_ = nullptr;
};

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_INSTRUMENT_H_
