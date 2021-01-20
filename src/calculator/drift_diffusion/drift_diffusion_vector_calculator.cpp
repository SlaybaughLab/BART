#include "calculator/drift_diffusion/drift_diffusion_vector_calculator.hpp"

#include "calculator/drift_diffusion/factory.hpp"

namespace bart::calculator::drift_diffusion {

template <int dim>
bool DriftDiffusionVectorCalculator<dim>::is_registered_ =
    DriftDiffusionVectorCalculatorIFactory<dim>::get()
        .RegisterConstructor(DriftDiffusionVectorCalculatorName::kDefaultImplementation,
                             []() -> std::unique_ptr<DriftDiffusionVectorCalculatorI<dim>>
                             { return std::make_unique<DriftDiffusionVectorCalculator<dim>>(); });

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
