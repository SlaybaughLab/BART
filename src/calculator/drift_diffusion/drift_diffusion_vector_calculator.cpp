#include "calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp"

namespace bart::calculator::drift_diffusion {

namespace  {
std::string error_preamble{ "Error in call to DriftDiffusion: " };
} // namespace


template<int dim>
auto DriftDiffusionVectorCalculator<dim>::DriftDiffusion(const double scalar_flux,
                                                         const double integrated_angular_flux,
                                                         const Tensor& shape_gradient,
                                                         const double sigma_t,
                                                         const double diffusion_coefficient) const -> Tensor {
  AssertThrow(integrated_angular_flux >= 0, dealii::ExcMessage(error_preamble + "integrated angular flux is negative"))
  AssertThrow(sigma_t > 0, dealii::ExcMessage(error_preamble + "possible called in void, sigma_t = 0"))
  if (scalar_flux == 0)
    return Tensor();
  return (integrated_angular_flux * shape_gradient / sigma_t
      - diffusion_coefficient * shape_gradient * scalar_flux) / scalar_flux;
}

template class DriftDiffusionVectorCalculator<1>;
template class DriftDiffusionVectorCalculator<2>;
template class DriftDiffusionVectorCalculator<3>;

} // namespace bart::calculator::drift_diffusion
