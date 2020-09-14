#include "calculator/fourier/fourier_transform_fftw.h"
#include "test_helpers/gmock_wrapper.h"

#include <cmath>
#include <complex.h>

namespace  {

namespace fftw {
#include <fftw3.h>
}

using namespace bart;

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

TEST_F(CalculatorFourierTransformFFTWTest, Cosine) {
  const double pi = M_PI;
  const int n = this->n_points;
  const double n_dob = n;
  //dealii::Vector<double> function(n);
  std::vector<std::complex<double>> function(n);
  for (int i = 0; i < function.size(); ++i) {
    const double x = 2*pi*static_cast<double>(i)/n_dob;
    function.at(i) = std::sin(x) * std::cos(2*pi*x);
  }

  std::vector<std::complex<double>> fourier_transform(n);

  fftw::fftw_complex* input = reinterpret_cast<fftw::fftw_complex*>(function.data());
  fftw::fftw_complex* output = reinterpret_cast<fftw::fftw_complex*>(fourier_transform.data());

  auto plan = fftw::fftw_plan_dft_1d(n, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw::fftw_execute(plan);
  fftw::fftw_destroy_plan(plan);

  std::vector<std::complex<double>> restored_function(n);
  fftw::fftw_complex* restored = reinterpret_cast<fftw::fftw_complex*>(restored_function.data());

  auto backwards_plan = fftw::fftw_plan_dft_1d(n, output, restored, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw::fftw_execute(backwards_plan);
  fftw::fftw_destroy_plan(backwards_plan);

  for (int i = 0; i < function.size(); ++i) {
    const auto& value = restored_function.at(i);
    EXPECT_NEAR(value.imag()/n_dob, 0, 1e-15);
    EXPECT_NEAR(value.real()/n_dob - function.at(i).real(), 0, 1e-15);
  }
}

} // namespace
