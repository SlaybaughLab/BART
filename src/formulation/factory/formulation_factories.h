#ifndef BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_
#define BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_

#include <memory>

#include "formulation/formulation_types.h"

// Formulation
#include "formulation/angular/self_adjoint_angular_flux_i.h"

// Dependencies
#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "quadrature/quadrature_set_i.h"

namespace bart {

namespace formulation {

namespace factory {

template <int dim>
std::unique_ptr<angular::SelfAdjointAngularFluxI<dim>> MakeSAAFFormulationPtr(
    const std::shared_ptr<domain::finite_element::FiniteElementI<dim>>&,
    const std::shared_ptr<data::CrossSections>&,
    const std::shared_ptr<quadrature::QuadratureSetI<dim>>&,
    const SAAFFormulationImpl implementation = SAAFFormulationImpl::kDefault);

} // namespace factory

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_FACTORY_FORMULATION_FACTORIES_H_
