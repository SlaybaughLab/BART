#ifndef BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_
#define BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_

#include <complex>
#include <memory>
#include <vector>

#include "calculator/fourier/fourier_transform_i.h"
#include "instrumentation/converter/converter_i.h"
#include "utility/has_dependencies.h"

namespace bart {

namespace instrumentation {

namespace converter {

namespace fourier {

using ComplexVector = std::vector<std::complex<double>>;

class FourierTransform : public ConverterI<ComplexVector, ComplexVector>,
                         public utility::HasDependencies {
 public:
  using FourierCalculator = calculator::fourier::FourierTransformI;

  explicit FourierTransform(std::unique_ptr<FourierCalculator>);
  ComplexVector Convert(const ComplexVector &input) const override;

  FourierCalculator* fourier_calculator_ptr() {
    return fourier_calculator_ptr_.get(); }
 private:
  std::unique_ptr<FourierCalculator> fourier_calculator_ptr_{nullptr};
};

} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_
