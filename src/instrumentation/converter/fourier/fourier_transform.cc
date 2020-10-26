#include "instrumentation/converter/fourier/fourier_transform.h"
#include "instrumentation/converter/factory.hpp"

namespace bart {

namespace instrumentation {

namespace converter {

namespace fourier {

FourierTransform::FourierTransform(
    std::unique_ptr<FourierCalculator> fourier_calculator_ptr)
    : FourierTransform(std::move(fourier_calculator_ptr), Normalized(true)) {}

FourierTransform::FourierTransform(
    std::unique_ptr<FourierCalculator> fourier_calculator_ptr,
    Normalized return_normalized)
    : fourier_calculator_ptr_(std::move(fourier_calculator_ptr)),
      returns_normalized_(return_normalized) {
  AssertPointerNotNull(fourier_calculator_ptr_.get(),
                       "Fourier transform calculator",
                       "Fourier transform converter constructor");
}

ComplexVector FourierTransform::Convert(const ComplexVector &input) const {
  return fourier_calculator_ptr_->CalculateDFT(input,
                                               Normalized(returns_normalized_));
}

bool FourierTransform::is_registered_ =
    ConverterIFactory<ComplexVector, ComplexVector, std::unique_ptr<FourierCalculator>, Normalized>::get()
    .RegisterConstructor(
        ConverterName::kFourierTransform,
        [](std::unique_ptr<FourierCalculator> fourier_calculator_ptr,
           Normalized return_normalized) {
          std::unique_ptr<ConverterI<ComplexVector, ComplexVector>> return_ptr =
              std::make_unique<FourierTransform>(std::move(fourier_calculator_ptr),
                                                 return_normalized);
          return return_ptr;
        }
    );


} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart