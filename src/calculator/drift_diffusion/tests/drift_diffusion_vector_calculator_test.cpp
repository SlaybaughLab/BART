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
  const double scalar_flux_{ 2.0 };
  const double diffision_coefficient_{ 3.0 };
  Tensor shape_gradient_, current_, expected_result_;

  auto SetUp() -> void override;
};

template <typename DimensionWrapper>
auto DriftDiffusionVectorCalculatorTest<DimensionWrapper>::SetUp() -> void {
  for (int i = 0; i < dim; ++i) {
    shape_gradient_[i] = 10 * (i + 1);
    current_[i] = i + 1;
    expected_result_[i] = (i + 1) * 15.5;
  }
}

TYPED_TEST_SUITE(DriftDiffusionVectorCalculatorTest, bart::testing::AllDimensions);

TYPED_TEST(DriftDiffusionVectorCalculatorTest, DefaultParameters) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  const auto result = this->test_calculator_.DriftDiffusionVector(this->scalar_flux_,
                                                                  this->current_,
                                                                  this->shape_gradient_,
                                                                  this->diffision_coefficient_);
  EXPECT_TRUE(AreEqual(this->expected_result_, result));
}

TYPED_TEST(DriftDiffusionVectorCalculatorTest, ZeroScalarFlux) {
  constexpr int dim{ this->dim };
  using Tensor = dealii::Tensor<1, dim>;
  const double scalar_flux{ 0 };
const auto result = this->test_calculator_.DriftDiffusionVector(scalar_flux,
                                                                this->current_,
                                                                this->shape_gradient_,
                                                                this->diffision_coefficient_);
  Tensor expected_result;
  expected_result = 0;
  EXPECT_TRUE(AreEqual(expected_result, result));
}

} // namespace
