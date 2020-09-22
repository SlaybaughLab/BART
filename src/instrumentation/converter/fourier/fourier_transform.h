#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_

#include <complex>
#include <vector>

#include "instrumentation/converter/converter_i.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace fourier {

using ComplexVector = std::vector<std::complex<double>>;

class FourierTransform : public ConverterI<ComplexVector, ComplexVector> {
 public:
  ComplexVector Convert(const ComplexVector &input) const override {
    return bart::instrumentation::converter::fourier::ComplexVector();
  }
};

} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_
