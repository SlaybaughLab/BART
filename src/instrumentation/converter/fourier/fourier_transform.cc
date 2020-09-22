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
  return bart::instrumentation::converter::fourier::ComplexVector();
}

} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart