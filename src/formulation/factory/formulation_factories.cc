#include "formulation/factory/formulation_factories.h"

#include "formulation/stamper.hpp"
#include "formulation/angular/self_adjoint_angular_flux.h"
#include "formulation/scalar/diffusion.h"

namespace bart {

namespace formulation {

namespace factory {

template <int dim>
std::unique_ptr<formulation::scalar::DiffusionI<dim>> MakeDiffusionPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>& finite_element_ptr,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>& cross_sections_ptr,
    const DiffusionFormulationImpl implementation) {
  std::unique_ptr<formulation::scalar::DiffusionI<dim>> return_ptr = nullptr;

  if (implementation == formulation::DiffusionFormulationImpl::kDefault) {
    return_ptr = std::move(
        std::make_unique<scalar::Diffusion<dim>>(
            finite_element_ptr, cross_sections_ptr));
  }

  return return_ptr;
}

template<int dim>
std::unique_ptr<updater::DiffusionUpdater<dim>> MakeDiffusionUpdater(
    std::unique_ptr<scalar::DiffusionI<dim>> formulation_ptr,
    std::unique_ptr<StamperI<dim>> stamper_ptr) {
  std::unique_ptr<updater::DiffusionUpdater<dim>> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<updater::DiffusionUpdater<dim>>(
          std::move(formulation_ptr),
          std::move(stamper_ptr)));

  return return_ptr;
}

template <int dim>
std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> MakeSAAFFormulationPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>& finite_element_ptr,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>& cross_sections_ptr,
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

template <int dim>
std::unique_ptr<formulation::updater::SAAFUpdater<dim>> MakeSAAFUpdater(
    std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> formulation_ptr,
    std::unique_ptr<StamperI<dim>> stamper_ptr,
    const std::shared_ptr<quadrature::QuadratureSetI<dim>>& quadrature_set_ptr) {
  std::unique_ptr<formulation::updater::SAAFUpdater<dim>> return_ptr = nullptr;

  return_ptr = std::move(
      std::make_unique<formulation::updater::SAAFUpdater<dim>>(
          std::move(formulation_ptr),
          std::move(stamper_ptr),
          quadrature_set_ptr));

  return return_ptr;
}

template <int dim>
std::unique_ptr<formulation::StamperI<dim>> MakeStamperPtr(
    const std::shared_ptr<domain::DomainI<dim>>& definition_ptr,
    const StamperImpl implementation) {
  std::unique_ptr<formulation::StamperI<dim>> return_ptr = nullptr;

  if (implementation == formulation::StamperImpl::kDefault) {
    return_ptr = std::move(
        std::make_unique<Stamper<dim>>(definition_ptr));
  }

  return return_ptr;
}

template std::unique_ptr<formulation::scalar::DiffusionI<1>> MakeDiffusionPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<1>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const DiffusionFormulationImpl);
template std::unique_ptr<formulation::scalar::DiffusionI<2>> MakeDiffusionPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<2>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const DiffusionFormulationImpl);
template std::unique_ptr<formulation::scalar::DiffusionI<3>> MakeDiffusionPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<3>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const DiffusionFormulationImpl);

template std::unique_ptr<updater::DiffusionUpdater<1>> MakeDiffusionUpdater<1>(
    std::unique_ptr<scalar::DiffusionI<1>>,
    std::unique_ptr<StamperI<1>>);
template std::unique_ptr<updater::DiffusionUpdater<2>> MakeDiffusionUpdater<2>(
    std::unique_ptr<scalar::DiffusionI<2>>,
    std::unique_ptr<StamperI<2>>);
template std::unique_ptr<updater::DiffusionUpdater<3>> MakeDiffusionUpdater<3>(
    std::unique_ptr<scalar::DiffusionI<3>>,
    std::unique_ptr<StamperI<3>>);

template std::unique_ptr<angular::SelfAdjointAngularFluxI<1>> MakeSAAFFormulationPtr<1>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<1>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<1>>&,
    const SAAFFormulationImpl);
template std::unique_ptr<angular::SelfAdjointAngularFluxI<2>> MakeSAAFFormulationPtr<2>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<2>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<2>>&,
    const SAAFFormulationImpl);
template std::unique_ptr<angular::SelfAdjointAngularFluxI<3>> MakeSAAFFormulationPtr<3>(
    const std::shared_ptr<domain::finite_element::FiniteElementI<3>>&,
    const std::shared_ptr<data::cross_sections::MaterialCrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<3>>&,
    const SAAFFormulationImpl);


template std::unique_ptr<formulation::updater::SAAFUpdater<1>> MakeSAAFUpdater<1>(
    std::unique_ptr<angular::SelfAdjointAngularFluxI<1>>,
    std::unique_ptr<StamperI<1>>,
    const std::shared_ptr<quadrature::QuadratureSetI<1>>&);
template std::unique_ptr<formulation::updater::SAAFUpdater<2>> MakeSAAFUpdater<2>(
    std::unique_ptr<angular::SelfAdjointAngularFluxI<2>>,
    std::unique_ptr<StamperI<2>>,
    const std::shared_ptr<quadrature::QuadratureSetI<2>>&);
template std::unique_ptr<formulation::updater::SAAFUpdater<3>> MakeSAAFUpdater<3>(
    std::unique_ptr<angular::SelfAdjointAngularFluxI<3>>,
    std::unique_ptr<StamperI<3>>,
    const std::shared_ptr<quadrature::QuadratureSetI<3>>&);

template std::unique_ptr<formulation::StamperI<1>> MakeStamperPtr(
    const std::shared_ptr<domain::DomainI<1>>&,
    const StamperImpl implementation);
template std::unique_ptr<formulation::StamperI<2>> MakeStamperPtr(
    const std::shared_ptr<domain::DomainI<2>>&,
    const StamperImpl implementation);
template std::unique_ptr<formulation::StamperI<3>> MakeStamperPtr(
    const std::shared_ptr<domain::DomainI<3>>&,
    const StamperImpl implementation);

} // namespace factory

} // namespace formulation

} // namespace bart
