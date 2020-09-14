#include "calculator/fourier/fourier_transform_fftw.h"
#include "test_helpers/gmock_wrapper.h"

#include <cmath>
#include <complex.h>

namespace  {

using namespace bart;
using ::testing::ContainerEq;
namespace fftw = bart::calculator::fourier::fftw;

class CalculatorFourierTransformFFTWTest : public ::testing::Test {
 public:
  using FourierTransformType = calculator::fourier::FourierTransformFFTW;
  void SetUp() override;
  static constexpr int n_points{1000};
  std::vector<std::complex<double>> function_, expected_fourier_transform_,
      normalized_expected_fourier_transform_;
  std::unique_ptr<FourierTransformType> test_transformer_ptr_;
};

void CalculatorFourierTransformFFTWTest::SetUp() {
  test_transformer_ptr_ = std::make_unique<FourierTransformType>(n_points);
  function_.resize(n_points);
  expected_fourier_transform_.resize(n_points);
  normalized_expected_fourier_transform_.resize(n_points);

  auto forward_input = reinterpret_cast<fftw::fftw_complex*>(function_.data());
  auto forward_output = reinterpret_cast<fftw::fftw_complex*>(expected_fourier_transform_.data());
  auto forwards_plan = fftw::fftw_plan_dft_1d(n_points, forward_input,
                                              forward_output, FFTW_FORWARD,
                                              FFTW_ESTIMATE);
  for (int i = 0; i < n_points; ++i) {
    const double pi = M_PI;
    const double x = 2*pi*static_cast<double>(i)/n_points;
    function_.at(i) = std::sin(x) * std::cos(2*pi*x);
  }
  fftw::fftw_execute(forwards_plan);
  fftw::fftw_destroy_plan(forwards_plan);

  normalized_expected_fourier_transform_ = expected_fourier_transform_;
  for (auto& value : normalized_expected_fourier_transform_)
    value /= n_points;
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

  auto calculated_fourier_transform =
      test_transformer_ptr_->CalculateDFT(function_);

  EXPECT_THAT(calculated_fourier_transform,
              ContainerEq(expected_fourier_transform_));
}

TEST_F(CalculatorFourierTransformFFTWTest, CosineStdVectorNormalized) {
  using bart::calculator::fourier::Normalized;

  auto calculated_fourier_transform =
      test_transformer_ptr_->CalculateDFT(function_, Normalized(true));

  EXPECT_THAT(calculated_fourier_transform,
              ContainerEq(normalized_expected_fourier_transform_));
}



} // namespace
