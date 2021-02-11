#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_HPP_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_HPP_

#include <complex>
#include <vector>

#include "calculator/fourier/fourier_transform_i.hpp"
#include "utility/named_type.h"

namespace bart::calculator::fourier {

namespace fftw {
#include <fftw3.h>
} // namespace fftw

class FourierTransformFFTW : public FourierTransformI {
 public:
  using FourierTransformI::ComplexVector;
  explicit FourierTransformFFTW(const int n_samples);
  ~FourierTransformFFTW();

  [[nodiscard]] auto CalculateDFT(const ComplexVector& input,
                                  Normalized normalized = Normalized(false)) -> ComplexVector override;
  [[nodiscard]] auto CalculateDFT(const dealii::Vector<double> &input,
                                  Normalized normalized = Normalized(false)) -> ComplexVector override;
  [[nodiscard]] auto n_samples() const -> int { return n_samples_; }
 private:
  const int n_samples_;
  std::vector<std::complex<double>> input_, output_;
  fftw::fftw_complex* input_ptr_;
  fftw::fftw_complex* output_ptr_;
  fftw::fftw_plan_s* plan_;
};

} // namespace bart::calculator::fourier

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_HPP_
