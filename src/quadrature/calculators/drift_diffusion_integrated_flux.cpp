#include "quadrature/calculators/drift_diffusion_integrated_flux.hpp"

namespace bart::quadrature::calculators {

template<int dim>
DriftDiffusionIntegratedFlux<dim>::DriftDiffusionIntegratedFlux(std::shared_ptr<QuadratureSet> quadrature_set_ptr)
    : quadrature_set_ptr_(quadrature_set_ptr) {
  this->AssertPointerNotNull(quadrature_set_ptr_.get(), "quadrature set", "DriftDiffusionIntegratedFlux constructor");
}

template class DriftDiffusionIntegratedFlux<1>;
template class DriftDiffusionIntegratedFlux<2>;
template class DriftDiffusionIntegratedFlux<3>;

} // namespace bart::quadrature::calculators
