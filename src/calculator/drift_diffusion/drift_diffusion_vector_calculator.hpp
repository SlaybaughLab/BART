#ifndef BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_HPP_
#define BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_HPP_

#include "calculator/drift_diffusion/drift_diffusion_vector_calculator_i.hpp"

namespace bart::calculator::drift_diffusion {

template <int dim>
class DriftDiffusionVectorCalculator : public DriftDiffusionVectorCalculatorI<dim> {
 public:
  using typename DriftDiffusionVectorCalculatorI<dim>::Tensor;

  [[nodiscard]] auto DriftDiffusionVector(const double scalar_flux,
                                          const Tensor& current,
                                          const Tensor& shape_gradient,
                                          const double diffusion_coefficient) const -> Tensor override;
 private:
  static bool is_registered_;
};

} // namespace bart::calculator::drift_diffusion

#endif //BART_SRC_CALCULATOR_DRIFT_DIFFUSION_DRIFT_DIFFUSION_VECTOR_CALCULATOR_HPP_
