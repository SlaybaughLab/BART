#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_HPP_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_HPP_

#include <complex>
#include <vector>

#include "calculator/fourier/fourier_transform_i.hpp"
#include "utility/named_type.h"

namespace bart::calculator::fourier {
//! Wrapper namespace for the C FFTW library.
namespace fftw {
#include <fftw3.h>
} // namespace fftw

/*! \brief Default implementation of the calculation of the DFT using the FFTW library.
 *
 * Due to the way FFTW works, the size of the input vector must be known when the class is instantiated. This has the
 * downside of requiring _a priori_ knowledge of the vector size. The benefit is that this enables FFTW to identify an
 * optimized algorithm for calculating the DFT.
 *
 * */
class FourierTransformFFTW : public FourierTransformI {
 public:
  using FourierTransformI::ComplexVector;
  /*! \brief Constructor.
   *
   * @param n_samples size of the vector that will be provided as input to this class.
   */
  explicit FourierTransformFFTW(const int n_samples);
  ~FourierTransformFFTW();

  [[nodiscard]] auto CalculateDFT(const ComplexVector& input,
                                  Normalized normalized = Normalized(false)) -> ComplexVector override;
  [[nodiscard]] auto CalculateDFT(const dealii::Vector<double> &input,
                                  Normalized normalized = Normalized(false)) -> ComplexVector override;
  /*! \brief Get number of samples identified during construction. */
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
