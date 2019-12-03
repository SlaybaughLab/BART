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

template<int dim>
auto CFEMSelfAdjointAngularFlux<dim>::Initialize(
    const formulation::CellPtr<dim> &cell_ptr) -> InitializationToken {
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage("Error in CFEMSelfAdjointAngularFlux Initialize, "
                                 "cell pointer is invalid."))

  finite_element_ptr_->SetCell(cell_ptr);
  shape_squared_ = {};

  for (int quad_index = 0; quad_index < cell_quadrature_points_; ++quad_index) {
    formulation::FullMatrix shape_squared(cell_degrees_of_freedom_,
                                          cell_degrees_of_freedom_);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) =
            finite_element_ptr_->ShapeValue(i, quad_index) *
            finite_element_ptr_->ShapeValue(j, quad_index);
      }
    }
    shape_squared_.insert_or_assign(quad_index, shape_squared);
  }

  return InitializationToken();
}

template class CFEMSelfAdjointAngularFlux<1>;
template class CFEMSelfAdjointAngularFlux<2>;
template class CFEMSelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
