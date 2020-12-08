#ifndef BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_HPP_

#include "quadrature/calculators/drift_diffusion_integrated_flux_i.hpp"

#include "quadrature/quadrature_set_i.h"
#include "utility/has_dependencies.h"

namespace bart::quadrature::calculators {

template <int dim>
class DriftDiffusionIntegratedFlux : public DriftDiffusionIntegratedFluxI, public utility::HasDependencies {
 public:
  using QuadratureSet = typename quadrature::QuadratureSetI<dim>;

  DriftDiffusionIntegratedFlux(std::shared_ptr<QuadratureSet>);

  auto quadrature_set_ptr() const -> QuadratureSet* { return quadrature_set_ptr_.get(); }

 private:
  std::shared_ptr<QuadratureSet> quadrature_set_ptr_{ nullptr };
};

} // namespace bart::quadrature::calculators

#endif //BART_SRC_QUADRATURE_CALCULATORS_DRIFT_DIFFUSION_INTEGRATED_FLUX_HPP_
