#include "formulation/scalar/drift_diffusion.hpp"

namespace bart::formulation::scalar {

template<int dim>
DriftDiffusion<dim>::DriftDiffusion(std::shared_ptr<FiniteElement> finite_element_ptr,
                                    std::shared_ptr<CrossSections> cross_sections_ptr,
                                    std::shared_ptr<DriftDiffusionCalculator> drift_diffusion_calculator_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      drift_diffusion_calculator_ptr_(drift_diffusion_calculator_ptr) {
  std::string function_name{"formulation::scalar::DriftDiffusion constructor"};
  AssertPointerNotNull(finite_element_ptr_.get(), "finite_element_ptr", function_name);
  AssertPointerNotNull(cross_sections_ptr_.get(), "cross_sections_ptr", function_name);
  AssertPointerNotNull(drift_diffusion_calculator_ptr_.get(), "drift_diffusion_calculator_ptr", function_name);
}

template class DriftDiffusion<1>;
template class DriftDiffusion<2>;
template class DriftDiffusion<3>;

} // namespace bart::formulation::scalar
