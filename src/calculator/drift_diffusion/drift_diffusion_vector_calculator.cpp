#include "calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp"

namespace bart::calculator::drift_diffusion {

template<int dim>
auto DriftDiffusionVectorCalculator<dim>::DriftDiffusionVector(const double scalar_flux,
                                                               const Tensor& current,
                                                               const Tensor& shape_gradient,
                                                               const double diffusion_coefficient) const -> Tensor {
  if (scalar_flux == 0) return Tensor();
  return (current + diffusion_coefficient * shape_gradient) / scalar_flux;
}

template class DriftDiffusionVectorCalculator<1>;
template class DriftDiffusionVectorCalculator<2>;
template class DriftDiffusionVectorCalculator<3>;

} // namespace bart::calculator::drift_diffusion
