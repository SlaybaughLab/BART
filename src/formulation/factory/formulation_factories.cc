#include "formulation/factory/formulation_factories.h"

#include "formulation/angular/self_adjoint_angular_flux.h"

namespace bart {

namespace formulation {

namespace factory {

template <int dim>
std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> MakeSAAFFormulationPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>& finite_element_ptr,
    const std::shared_ptr<data::CrossSections>& cross_sections_ptr,
    const std::shared_ptr<quadrature::QuadratureSetI<dim>>& quadrature_set_ptr,
    const SAAFFormulationImpl implementation) {
  std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> return_ptr = nullptr;

  if (implementation == formulation::SAAFFormulationImpl::kDefault) {
    return_ptr = std::move(
        std::make_unique<angular::SelfAdjointAngularFlux<dim>>(
            finite_element_ptr, cross_sections_ptr, quadrature_set_ptr));
  }

  return return_ptr;
}

template std::unique_ptr<angular::SelfAdjointAngularFluxI<1>> MakeSAAFFormulationPtr<1>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<1>>&,
    const std::shared_ptr<data::CrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<1>>&,
    const SAAFFormulationImpl);
template std::unique_ptr<angular::SelfAdjointAngularFluxI<2>> MakeSAAFFormulationPtr<2>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<2>>&,
    const std::shared_ptr<data::CrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<2>>&,
    const SAAFFormulationImpl);
template std::unique_ptr<angular::SelfAdjointAngularFluxI<3>> MakeSAAFFormulationPtr<3>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<3>>&,
    const std::shared_ptr<data::CrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<3>>&,
    const SAAFFormulationImpl);
} // namespace factory

} // namespace formulation

} // namespace bart
