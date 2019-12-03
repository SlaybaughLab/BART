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
  omega_dot_gradient_ = {};

  /* Precalculated values are held in maps that are indexed by the cell
   * quadrature point where they are valid */
  for (int cell_quad_index = 0; cell_quad_index < cell_quadrature_points_;
       ++cell_quad_index) {
    formulation::FullMatrix shape_squared(cell_degrees_of_freedom_,
                                          cell_degrees_of_freedom_);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) =
            finite_element_ptr_->ShapeValue(i, cell_quad_index) *
            finite_element_ptr_->ShapeValue(j, cell_quad_index);
      }
    }
    shape_squared_.insert_or_assign(cell_quad_index, shape_squared);

    /* Each cell quadrature point has a further mapping based on the angular
     * quadrature point. */
    std::map<int, FullMatrix> omega_dot_gradient_squared_map;

    for (int angle_index : quadrature_set_ptr_->quadrature_point_indices()) {
      std::vector<double> angle_omega_dot_gradient(cell_degrees_of_freedom_);
      dealii::Vector<double> angle_omega_dot_gradient_vector(
          cell_degrees_of_freedom_);

      auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(
          quadrature::QuadraturePointIndex(angle_index));
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        double entry_value = quadrature_point_ptr->cartesian_position_tensor() *
            finite_element_ptr_->ShapeGradient(i, cell_quad_index);
        angle_omega_dot_gradient.at(i) = entry_value;
        angle_omega_dot_gradient_vector[i] = entry_value;
      }

      omega_dot_gradient_.insert_or_assign({cell_quad_index, angle_index},
                                           angle_omega_dot_gradient);

      FullMatrix omega_dot_gradient_squared(cell_degrees_of_freedom_,
                                            cell_degrees_of_freedom_);
      omega_dot_gradient_squared.outer_product(angle_omega_dot_gradient_vector,
                                               angle_omega_dot_gradient_vector);
      omega_dot_gradient_squared_map.insert_or_assign(angle_index,
                                                      omega_dot_gradient_squared);
    }

    omega_dot_gradient_squared_.insert_or_assign(cell_quad_index,
                                                 omega_dot_gradient_squared_map);
  }


  return InitializationToken();
}

template class CFEMSelfAdjointAngularFlux<1>;
template class CFEMSelfAdjointAngularFlux<2>;
template class CFEMSelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
