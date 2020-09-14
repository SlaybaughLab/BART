#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_

#include <complex>
#include <vector>

#include "calculator/fourier/fourier_transform_i.h"
#include "utility/named_type.h"

namespace bart {

namespace calculator {

namespace fourier {

namespace fftw {
#include <fftw3.h>
} // namespace fftw

class FourierTransformFFTW : public FourierTransformI {
 public:

  explicit FourierTransformFFTW(const int n_samples);
  ~FourierTransformFFTW();

  std::vector<std::complex<double>> CalculateDFT(
      const std::vector<std::complex<double>>& input,
      Normalized normalized = Normalized(false)) override;
  std::vector<std::complex<double>> CalculateDFT(
      const dealii::Vector<double> &input,
      Normalized normalized = Normalized(false)) override;

  int n_samples() const { return n_samples_; }
 private:
  const int n_samples_;
  std::vector<std::complex<double>> input_, output_;
  fftw::fftw_complex* input_ptr_;
  fftw::fftw_complex* output_ptr_;
  fftw::fftw_plan_s* plan_;
};

} // namespace fourier

} // namespace calculator

} // namespace bart

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_
