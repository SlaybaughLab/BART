#include "calculator/fourier/fourier_transform_fftw.h"

namespace bart {

namespace calculator {

namespace fourier {

FourierTransformFFTW::FourierTransformFFTW(const int n_samples)
    : n_samples_(n_samples) {
  input_.resize(n_samples);
  output_.resize(n_samples);
  input_ptr_ = reinterpret_cast<fftw::fftw_complex*>(input_.data());
  output_ptr_ = reinterpret_cast<fftw::fftw_complex*>(output_.data());
  plan_ = fftw::fftw_plan_dft_1d(n_samples_, input_ptr_, output_ptr_,
                                 FFTW_FORWARD, FFTW_ESTIMATE_PATIENT);
}

FourierTransformFFTW::~FourierTransformFFTW() {
  fftw::fftw_destroy_plan(plan_);
}

std::vector<std::complex<double>> FourierTransformFFTW::CalculateDFT(
    const std::vector<std::complex<double>>& input,
    Normalized normalized) {
  input_ = input;
  fftw::fftw_execute(plan_);
  if (normalized) {
    for (auto& value : output_)
      value /= n_samples_;
  }
  return output_;
}
std::vector<std::complex<double>> FourierTransformFFTW::CalculateDFT(
    const dealii::Vector<double>& input,
    Normalized normalized) {
  std::vector<std::complex<double>> input_copy(n_samples_);
  for (int i = 0; i < n_samples_; ++i) {
    input_copy.at(i).real(input[i]);
  }
  return CalculateDFT(input_copy, normalized);
}

} // namespace fourier

} // namespace calculator

} // namespace bart
