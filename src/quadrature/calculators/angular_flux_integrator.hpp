#ifndef BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_HPP_
#define BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_HPP_

#include "quadrature/calculators/angular_flux_integrator_i.hpp"

#include "quadrature/quadrature_set_i.h"
#include "utility/has_dependencies.h"

namespace bart::quadrature::calculators {

template <int dim>
class AngularFluxIntegrator : public AngularFluxIntegratorI, public utility::HasDependencies {
 public:
  using QuadratureSet = typename quadrature::QuadratureSetI<dim>;

  AngularFluxIntegrator(std::shared_ptr<QuadratureSet>);

  [[nodiscard]] auto Integrate(const VectorMap&) const -> Vector override;

  auto quadrature_set_ptr() const -> QuadratureSet* { return quadrature_set_ptr_.get(); }

 private:
  std::shared_ptr<QuadratureSet> quadrature_set_ptr_{ nullptr };
};

} // namespace bart::quadrature::calculators

#endif //BART_SRC_QUADRATURE_CALCULATORS_ANGULAR_FLUX_INTEGRATOR_HPP_
