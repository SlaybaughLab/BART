#include "test_helpers/gmock_wrapper.h"

#include "calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp"
#include "test_helpers/test_assertions.hpp"

namespace  {

using namespace bart;
using test_helpers::AreEqual;

template <typename DimensionWrapper>
class DriftDiffusionVectorCalculatorTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using Tensor = typename dealii::Tensor<1, dim>;
  using TestCalculator = typename calculator::drift_diffusion::DriftDiffusionVectorCalculator<dim>;

  // Test object
  TestCalculator test_calculator_;

  // Test parameters
  const double scalar_flux_{ 5.0 };
  const double integrated_angular_flux_{ 10.2 };
  const double sigma_t_{ 0.25 };
  const double diffision_coefficient_{ 2.3 };
  Tensor shape_gradient;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto DriftDiffusionVectorCalculatorTest<DimensionWrapper>::SetUp() -> void {
  shape_gradient[0] = 1.1;
  if (dim > 1) {
    shape_gradient[1] = 2.2;
    if (dim > 2) {
      shape_gradient[2] = 3.3;
    }
  }
}

TYPED_TEST_SUITE(DriftDiffusionVectorCalculatorTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionVectorCalculatorTest, DefaultParameters) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  const auto result = this->test_calculator_.DriftDiffusion(this->scalar_flux_, this->integrated_angular_flux_,
                                                            this->shape_gradient, this->sigma_t_,
                                                            this->diffision_coefficient_);
  Tensor expected_result;
  expected_result[0] = 6.446;
  if (dim > 1) {
    expected_result[1] = 12.892;
    if (dim > 2) {
      expected_result[2] = 19.338;
    }
  }
  EXPECT_TRUE(AreEqual(expected_result, result));
}

TYPED_TEST(DriftDiffusionVectorCalculatorTest, ZeroSigmaT) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  const double sigma_t{ 0 };
  EXPECT_ANY_THROW({
  [[maybe_unused]] auto result = this->test_calculator_.DriftDiffusion(this->scalar_flux_,
                                                                       this->integrated_angular_flux_,
                                                                       this->shape_gradient,
                                                                       sigma_t,
                                                                       this->diffision_coefficient_);
                   });
}

TYPED_TEST(DriftDiffusionVectorCalculatorTest, NegativeAngularFlux) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  EXPECT_ANY_THROW({
    [[maybe_unused]] auto result = this->test_calculator_.DriftDiffusion(this->scalar_flux_,
                                                                         -this->integrated_angular_flux_,
                                                                         this->shape_gradient,
                                                                         this->sigma_t_,
                                                                         this->diffision_coefficient_);
                   });
}

TYPED_TEST(DriftDiffusionVectorCalculatorTest, ZeroScalarFlux) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  const double scalar_flux{ 0 };
  const auto result = this->test_calculator_.DriftDiffusion(scalar_flux, this->integrated_angular_flux_,
                                                            this->shape_gradient, this->sigma_t_,
                                                            this->diffision_coefficient_);
  Tensor expected_result;
  expected_result[0] = 0;
  if (dim > 1) {
    expected_result[1] = 0;
    if (dim > 2) {
      expected_result[2] = 0;
    }
  }
  EXPECT_TRUE(AreEqual(expected_result, result));
}

} // namespace
