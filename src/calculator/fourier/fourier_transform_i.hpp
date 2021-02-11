#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

#include "utility/named_type.h"

namespace bart {

namespace calculator {

namespace fourier {

using Normalized = utility::NamedType<bool, struct NormalizedStruct>;

class FourierTransformI {
 public:
  virtual ~FourierTransformI() = default;
  virtual std::vector<std::complex<double>> CalculateDFT(
      const std::vector<std::complex<double>>& input,
      Normalized) = 0;
  virtual std::vector<std::complex<double>> CalculateDFT(
      const dealii::Vector<double>& input,
      Normalized) = 0;
};

} // namespace fourier

} // namespace calculator

} // namespace bart

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_HPP_
