#ifndef BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_
#define BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_

#include <memory>

#include "formulation/formulation_types.hpp"

// Formulations
#include "formulation/angular/self_adjoint_angular_flux_i.h"
#include "formulation/scalar/diffusion_i.hpp"

// Stamper
#include "formulation/stamper_i.hpp"

// Updaters
#include "formulation/updater/saaf_updater.h"
#include "formulation/updater/diffusion_updater.hpp"

// Dependencies
#include "data/cross_sections/material_cross_sections.hpp"
#include "domain/domain_i.hpp"
#include "domain/finite_element/finite_element_i.hpp"
#include "quadrature/quadrature_set_i.hpp"

namespace bart {

namespace formulation {

namespace factory {

template <int dim>
std::unique_ptr<formulation::scalar::DiffusionI<dim>> MakeDiffusionPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>&,
    const std::shared_ptr<data::cross_sections::CrossSectionsI>&,
    const DiffusionFormulationImpl implementation = DiffusionFormulationImpl::kDefault);

template <int dim>
std::unique_ptr<updater::DiffusionUpdater<dim>> MakeDiffusionUpdater(
    std::unique_ptr<scalar::DiffusionI<dim>>,
    std::unique_ptr<StamperI<dim>>);

template <int dim>
std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> MakeSAAFFormulationPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>&,
    const std::shared_ptr<data::cross_sections::CrossSectionsI>&,
    const std::shared_ptr<quadrature::QuadratureSetI<dim>>&,
    const SAAFFormulationImpl implementation = SAAFFormulationImpl::kDefault);

template <int dim>
std::unique_ptr<formulation::updater::SAAFUpdater<dim>> MakeSAAFUpdater(
    std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>>,
    std::unique_ptr<StamperI<dim>>,
    const std::shared_ptr<quadrature::QuadratureSetI<dim>>&);

template <int dim>
std::unique_ptr<formulation::StamperI<dim>> MakeStamperPtr(
    const std::shared_ptr<domain::DomainI<dim>>&,
    const StamperImpl implementation = StamperImpl::kDefault);


} // namespace factory

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_
