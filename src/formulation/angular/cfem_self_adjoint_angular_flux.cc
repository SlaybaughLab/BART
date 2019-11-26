#include "formulation/angular/cfem_self_adjoint_angular_flux.h"

namespace bart {

namespace formulation {

namespace angular {

template <int dim>
CFEMSelfAdjointAngularFlux<dim>::CFEMSelfAdjointAngularFlux(
    std::shared_ptr<domain::finite_element::FiniteElementI<dim>> finite_element_ptr,
    std::shared_ptr<data::CrossSections> cross_sections_ptr,
    std::shared_ptr<quadrature::QuadratureSetI<dim>> quadrature_set_ptr)
    : finite_element_ptr_(finite_element_ptr),
      cross_sections_ptr_(cross_sections_ptr),
      quadrature_set_ptr_(quadrature_set_ptr),
      cell_degrees_of_freedom_(finite_element_ptr->dofs_per_cell()),
      cell_quadrature_points_(finite_element_ptr->n_cell_quad_pts()),
      face_quadrature_points_(finite_element_ptr->n_face_quad_pts()) {}

template class CFEMSelfAdjointAngularFlux<1>;
template class CFEMSelfAdjointAngularFlux<2>;
template class CFEMSelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
