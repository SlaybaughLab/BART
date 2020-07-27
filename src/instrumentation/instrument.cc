#include "instrumentation/instrument.h"

#include <deal.II/base/exceptions.h>
#include <string>

#include "convergence/status.h"

namespace bart {

namespace instrumentation {

template<typename InputType, typename OutputType>
Instrument<InputType, OutputType>::Instrument(
    std::unique_ptr<ConverterType> converter_ptr,
    std::unique_ptr<OutputterType> outputter_ptr)
    : converter_ptr_(std::move(converter_ptr)),
      outputter_ptr_(std::move(outputter_ptr)) {
  AssertThrow(converter_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of instrument, "
                                 "converter_ptr passed is null"))
  AssertThrow(outputter_ptr_ != nullptr,
              dealii::ExcMessage("Error in constructor of instrument, "
                                 "outputter_ptr passed is null"))
}

template<typename InputType, typename OutputType>
void Instrument<InputType, OutputType>::Read(const InputType &input) {
  outputter_ptr_->Output(converter_ptr_->Convert(input));
}

template class Instrument<convergence::Status, std::string>;

} // namespace instrumentation

} // namespace bart
