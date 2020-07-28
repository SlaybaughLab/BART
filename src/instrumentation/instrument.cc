#include "instrumentation/instrument.h"

#include <deal.II/base/exceptions.h>
#include <string>

#include "convergence/status.h"

namespace bart {

namespace instrumentation {

template <typename InputType, typename OutputType>
Instrument<InputType, OutputType>::Instrument(
    std::unique_ptr<ConverterType> converter_ptr,
    std::unique_ptr<OutputterType> outputter_ptr)
    : converter_ptr_(std::move(converter_ptr)),
      outputter_ptr_(std::move(outputter_ptr)) {
  AssertNotNull(converter_ptr_.get(), "converter", __FUNCTION__ );
  AssertNotNull(outputter_ptr_.get(), "outputter", __FUNCTION__ );
}

template<typename InputType, typename OutputType>
void Instrument<InputType, OutputType>::Read(const InputType &input) {
  outputter_ptr_->Output(converter_ptr_->Convert(input));
}

template<typename InputType, typename OutputType>
template<typename T>
void Instrument<InputType, OutputType>::AssertNotNull(T* ptr,
                                                      std::string dependency,
                                                      std::string function_name) {
  AssertThrow(ptr != nullptr,
              dealii::ExcMessage("Error in " + function_name + " passed "
                                     + dependency + " pointer is null"))
}

template class Instrument<convergence::Status, std::string>;

} // namespace instrumentation

} // namespace bart
