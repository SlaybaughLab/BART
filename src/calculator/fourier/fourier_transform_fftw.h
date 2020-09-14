#ifndef BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_
#define BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_

#include "calculator/fourier/fourier_transform_i.h"

namespace bart {

namespace calculator {

namespace fourier {

class FourierTransformFFTW : public FourierTransformI {
 public:
  FourierTransformFFTW() = default;
};

} // namespace fourier

} // namespace calculator

} // namespace bart

#endif //BART_SRC_CALCULATOR_FOURIER_FOURIER_TRANSFORM_FFTW_H_
