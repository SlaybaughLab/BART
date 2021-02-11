#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

#include "utility/named_type.h"

namespace bart::calculator::fourier {

using Normalized = utility::NamedType<bool, struct NormalizedStruct>;

class FourierTransformI {
 public:
  using ComplexVector = std::vector<std::complex<double>>;
  virtual ~FourierTransformI() = default;
  virtual auto CalculateDFT(const ComplexVector& input, Normalized) -> ComplexVector = 0;
  virtual auto CalculateDFT(const dealii::Vector<double>& input, Normalized) -> ComplexVector = 0;
};

} // namespace bart::calculator::fourier

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_
