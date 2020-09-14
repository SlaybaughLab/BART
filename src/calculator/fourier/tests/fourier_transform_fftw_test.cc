#include "calculator/fourier/fourier_transform_fftw.h"
#include "test_helpers/gmock_wrapper.h"

#include <cmath>
#include <complex.h>

namespace  {

using namespace bart;
namespace fftw = bart::calculator::fourier::fftw;

class CalculatorFourierTransformFFTWTest : public ::testing::Test {
 public:
  using FourierTransformType = calculator::fourier::FourierTransformFFTW;
  void SetUp() override;
  static constexpr int n_points{1000};
  std::unique_ptr<FourierTransformType> test_transformer_ptr_;
};

void CalculatorFourierTransformFFTWTest::SetUp() {
  test_transformer_ptr_ = std::make_unique<FourierTransformType>(n_points);
}

TEST_F(CalculatorFourierTransformFFTWTest, Constructor) {
  EXPECT_NO_THROW({
    FourierTransformType test_transformer(n_points);
                  });
}

TEST_F(CalculatorFourierTransformFFTWTest, Getters) {
  EXPECT_EQ(test_transformer_ptr_->n_samples(), n_points);
}

TEST_F(CalculatorFourierTransformFFTWTest, CosineStdVector) {
  const double pi = M_PI;
  const int n = this->n_points;
  const double n_dob = n;
  std::vector<std::complex<double>> function(n);
  std::vector<std::complex<double>> expected_fourier_transform(n);

  fftw::fftw_complex* forward_input = reinterpret_cast<fftw::fftw_complex*>(function.data());
  fftw::fftw_complex* forward_output = reinterpret_cast<fftw::fftw_complex*>(expected_fourier_transform.data());

  auto forwards_plan = fftw::fftw_plan_dft_1d(n, forward_input, forward_output, FFTW_FORWARD, FFTW_ESTIMATE);

  for (int i = 0; i < n; ++i) {
    const double x = 2*pi*static_cast<double>(i)/n_dob;
    function.at(i) = std::sin(x) * std::cos(2*pi*x);
  }

  //expected_fourier_transform = test_transformer_ptr_->CalculateDFT(function);
  auto calculated_fourier_transform = test_transformer_ptr_->CalculateDFT(function);
  fftw::fftw_execute(forwards_plan);
  fftw::fftw_destroy_plan(forwards_plan);

  for (int i = 0; i < n; ++i) {
    EXPECT_EQ(calculated_fourier_transform.at(i),
              expected_fourier_transform.at(i));
  }

}

} // namespace
