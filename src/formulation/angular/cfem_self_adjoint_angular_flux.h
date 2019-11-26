#ifndef BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
#define BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_

#include "data/cross_sections.h"
#include "domain/finite_element/finite_element_i.h"
#include "formulation/angular/cfem_self_adjoint_angular_flux_i.h"
#include "quadrature/quadrature_set_i.h"

#include <memory>

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
class CFEMSelfAdjointAngularFlux : public CFEMSelfAdjointAngularFluxI<dim> {
 public:
  CFEMSelfAdjointAngularFlux(
      std::shared_ptr<domain::finite_element::FiniteElementI<dim>>,
      std::shared_ptr<data::CrossSections>,
      std::shared_ptr<quadrature::QuadratureSetI<dim>>);

  domain::finite_element::FiniteElementI<dim>* finite_element_ptr() const {
    return finite_element_ptr_.get(); }
  data::CrossSections* cross_sections_ptr() const {
    return cross_sections_ptr_.get(); }
  quadrature::QuadratureSetI<dim>* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get(); }

 protected:
  std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr_;
  std::shared_ptr<data::CrossSections> cross_sections_ptr_;
  std::shared_ptr<quadrature::QuadratureSetI<dim>> quadrature_set_ptr_;
  const int cell_degrees_of_freedom_ = 0;
  const int cell_quadrature_points_ = 0;
  const int face_quadrature_points_ = 0;

};

} // namespace angular

} // namespace formulation

} //namespace bart

#endif //BART_SRC_FORMULATION_ANGULAR_CFEM_SELF_ADJOINT_ANGULAR_FLUX_H_
