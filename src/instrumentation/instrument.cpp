#include "instrumentation/instrument.h"

#include <deal.II/base/exceptions.h>
#include <string>
#include <deal.II/lac/vector.h>

#include "convergence/status.hpp"
#include "utility/colors.hpp"
#include "system/moments/spherical_harmonic_i.h"

namespace bart::instrumentation {

template <typename InputType, typename OutputType>
Instrument<InputType, OutputType>::Instrument(
    std::unique_ptr<ConverterType> converter_ptr,
    std::unique_ptr<OutstreamType> outstream_ptr)
    : converter_ptr_(std::move(converter_ptr)),
      outstream_ptr_(std::move(outstream_ptr)) {
  AssertNotNull(converter_ptr_.get(), "converter", __FUNCTION__ );
  AssertNotNull(outstream_ptr_.get(), "outstream", __FUNCTION__ );
}

template<typename InputType, typename OutputType>
void Instrument<InputType, OutputType>::Read(const InputType &input) {
  outstream_ptr_->Output(converter_ptr_->Convert(input));
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
template class Instrument<std::pair<int, double>, std::string>;
template class Instrument<std::pair<std::string, utility::Color>, std::string>;
template class Instrument<dealii::Vector<double>, std::string>;
template class Instrument<system::moments::SphericalHarmonicI, std::string>;
template class Instrument<double, std::string>;

} // namespace bart::instrumentation
