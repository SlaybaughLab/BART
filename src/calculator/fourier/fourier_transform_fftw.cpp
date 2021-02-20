#include "calculator/fourier/fourier_transform_fftw.hpp"

namespace bart::calculator::fourier {

FourierTransformFFTW::FourierTransformFFTW(const int n_samples)
    : n_samples_(n_samples) {
  input_.resize(n_samples);
  output_.resize(n_samples);
  input_ptr_ = reinterpret_cast<fftw::fftw_complex*>(input_.data());
  output_ptr_ = reinterpret_cast<fftw::fftw_complex*>(output_.data());
  plan_ = fftw::fftw_plan_dft_1d(n_samples_, input_ptr_, output_ptr_, FFTW_FORWARD, FFTW_ESTIMATE_PATIENT);
  this->set_description("fourier transform calculator using FFTW", utility::DefaultImplementation(true));
}

FourierTransformFFTW::~FourierTransformFFTW() {
  fftw::fftw_destroy_plan(plan_);
}

auto FourierTransformFFTW::CalculateDFT(const ComplexVector& input, Normalized normalized) -> ComplexVector {
  input_ = input;
  fftw::fftw_execute(plan_);
  if (normalized) {
    for (auto& value : output_)
      value /= n_samples_;
  }
  return output_;
}

auto FourierTransformFFTW::CalculateDFT(const dealii::Vector<double>& input, Normalized normalized) -> ComplexVector {
  std::vector<std::complex<double>> input_copy(n_samples_);
  for (int i = 0; i < n_samples_; ++i)
    input_copy.at(i).real(input[i]);
  return CalculateDFT(input_copy, normalized);
}

} // namespace bart::calculator::fourier