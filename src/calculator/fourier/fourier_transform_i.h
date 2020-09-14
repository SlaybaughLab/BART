#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_

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
};

} // namespace fourier

} // namespace calculator

} // namespace bart

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_
