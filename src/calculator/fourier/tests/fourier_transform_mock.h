#ifndef BART_SRC_CALCULATOR_FOURIER_TESTS_FOURIER_TRANSFORM_MOCK_H_
#define BART_SRC_CALCULATOR_FOURIER_TESTS_FOURIER_TRANSFORM_MOCK_H_

#include "calculator/fourier/fourier_transform_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace calculator {

namespace fourier {

class FourierTransformMock : public FourierTransformI {
 public:
  MOCK_METHOD(std::vector<std::complex<double>>,
              CalculateDFT,
              (const std::vector<std::complex<double>>&, Normalized),
              (override));
  MOCK_METHOD(std::vector<std::complex<double>>,
              CalculateDFT,
              (const dealii::Vector<double>&, Normalized),
              (override));
};


} // namespace fourier

} // namespace calculator

} // namespace bart



#endif //BART_SRC_CALCULATOR_FOURIER_TESTS_FOURIER_TRANSFORM_MOCK_H_
