#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_

#include <complex>
#include <vector>

#include <deal.II/lac/vector.h>

namespace bart {

namespace calculator {

namespace fourier {

class FourierTransformI {
 public:
  virtual ~FourierTransformI() = default;
  virtual std::vector<std::complex<double>> CalculateDFT(
      const std::vector<std::complex<double>>& input) = 0;
};

} // namespace fourier

} // namespace calculator

} // namespace bart

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_I_H_
