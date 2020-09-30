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
  using FourierCalculator = bart::calculator::fourier::FourierTransformI;
  using Normalized = bart::calculator::fourier::Normalized;

  explicit FourierTransform(std::unique_ptr<FourierCalculator>);
  FourierTransform(std::unique_ptr<FourierCalculator>, Normalized);
  ComplexVector Convert(const ComplexVector &input) const override;

  FourierCalculator* fourier_calculator_ptr() {
    return fourier_calculator_ptr_.get(); }
  FourierTransform& returns_normalized(bool to_set) {
    returns_normalized_ = to_set;
    return *this;}
  bool returns_normalized() const { return returns_normalized_; }
 private:
  std::unique_ptr<FourierCalculator> fourier_calculator_ptr_{nullptr};
  bool returns_normalized_{true};
  static bool is_registered_;
};

} // namespace fourier

} // namespace converter

} // namespace instrumentation

} // namespace bart

#endif //BART_SRC_INSTRUMENTATION_CONVERTER_FOURIER_FOURIER_TRANSFORM_H_
