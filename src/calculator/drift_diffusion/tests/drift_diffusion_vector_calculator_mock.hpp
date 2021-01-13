#ifndef BART_SRC_CALCULATOR_DRIFT_DIFFUSION_TESTS_DRIFT_DIFFUSION_VECTOR_CALCULATOR_MOCK_HPP_
#define BART_SRC_CALCULATOR_DRIFT_DIFFUSION_TESTS_DRIFT_DIFFUSION_VECTOR_CALCULATOR_MOCK_HPP_

#include "calculator/drift_diffusion/drift_diffusion_vector_calculator_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::drift_diffusion {

template <int dim>
class DriftDiffusionVectorCalculatorMock : public DriftDiffusionVectorCalculatorI<dim> {
 public:
  using typename DriftDiffusionVectorCalculatorI<dim>::Tensor;
  MOCK_METHOD(Tensor, DriftDiffusionVector, (const double scalar_flux, const Tensor& current,
      const Tensor& shape_gradient, const double diffusion_coefficient), (const, override));
};

} // namespace bart::calculator::drift_diffusion

#endif //BART_SRC_CALCULATOR_DRIFT_DIFFUSION_TESTS_DRIFT_DIFFUSION_VECTOR_CALCULATOR_MOCK_HPP_
