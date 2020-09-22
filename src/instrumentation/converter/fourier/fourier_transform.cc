#include "instrumentation/converter/fourier/fourier_transform.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace fourier {

FourierTransform::FourierTransform(
    std::unique_ptr<FourierCalculator> fourier_calculator_ptr)
    : fourier_calculator_ptr_(std::move(fourier_calculator_ptr)) {
  AssertPointerNotNull(fourier_calculator_ptr_.get(),
                       "Fourier transform calculator",
                       "Fourier transform converter constructor");
}

ComplexVector FourierTransform::Convert(const ComplexVector &input) const {
  using calculator::fourier::Normalized;
  return fourier_calculator_ptr_->CalculateDFT(input, Normalized(true));
}

} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart