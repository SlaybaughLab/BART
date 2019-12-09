#include "formulation/angular/cfem_self_adjoint_angular_flux.h"

#include <sstream>

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
    for (int angle_index : quadrature_set_ptr_->quadrature_point_indices()) {
      auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(
          quadrature::QuadraturePointIndex(angle_index));

      dealii::Vector<double> omega_dot_gradient_vector(cell_degrees_of_freedom_);
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        omega_dot_gradient_vector[i] =
            quadrature_point_ptr->cartesian_position_tensor() *
            finite_element_ptr_->ShapeGradient(i, cell_quad_index);
      }

      omega_dot_gradient_.insert_or_assign({cell_quad_index, angle_index},
                                           omega_dot_gradient_vector);

      FullMatrix omega_dot_gradient_squared(cell_degrees_of_freedom_,
                                            cell_degrees_of_freedom_);
      omega_dot_gradient_squared.outer_product(omega_dot_gradient_vector,
                                               omega_dot_gradient_vector);
      omega_dot_gradient_squared_.insert_or_assign(
          {cell_quad_index, angle_index},
          omega_dot_gradient_squared);
    }
  }
  return InitializationToken();
}

template<int dim>
void CFEMSelfAdjointAngularFlux<dim>::FillCellCollisionTerm(
    FullMatrix &to_fill,
    const InitializationToken,
    const CellPtr<dim> &cell_ptr,
    const system::EnergyGroup group_number) {

  ValidateMatrixSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const double sigma_t =
      cross_sections_ptr_->sigma_t.at(material_id).at(group_number.get());

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += sigma_t * finite_element_ptr_->ShapeValue(i, q) *
            finite_element_ptr_->ShapeValue(j, q) * jacobian;
      }
    }
  }
}

template<int dim>
void CFEMSelfAdjointAngularFlux<dim>::FillCellFixedSourceTerm(
    Vector &to_fill,
    const InitializationToken,
    const CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number) {
  ValidateVectorSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const double inverse_sigma_t =
      cross_sections_ptr_->inverse_sigma_t.at(material_id).at(group_number.get());
  const double q_per_ster =
      cross_sections_ptr_->q_per_ster.at(material_id).at(group_number.get());
  const int angle_index = quadrature_set_ptr_->GetQuadraturePointIndex(
      quadrature_point);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    const std::vector<double> omega_dot_gradient = OmegaDotGradient(q,
        quadrature::QuadraturePointIndex(angle_index));
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      to_fill(i) += jacobian * q_per_ster *
          (finite_element_ptr_->ShapeValue(i, q)
              + omega_dot_gradient.at(i) * inverse_sigma_t);
    }
  }
}

template<int dim>
void CFEMSelfAdjointAngularFlux<dim>::FillCellScatteringSourceTerm(
    Vector &to_fill,
    const InitializationToken token,
    const CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number,
    const system::moments::MomentVector &in_group_moment,
    const system::moments::MomentsMap &group_moments) {


}

template<int dim>
void CFEMSelfAdjointAngularFlux<dim>::FillCellStreamingTerm(
    FullMatrix &to_fill,
    const InitializationToken,
    const CellPtr<dim> &cell_ptr,
    const std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point,
    const system::EnergyGroup group_number) {

  ValidateMatrixSizeAndSetCell(cell_ptr, to_fill, __FUNCTION__);

  const int material_id = cell_ptr->material_id();
  const double inverse_sigma_t =
      cross_sections_ptr_->inverse_sigma_t.at(material_id).at(group_number.get());
  const int angle_index = quadrature_set_ptr_->GetQuadraturePointIndex(
      quadrature_point);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    const double jacobian = finite_element_ptr_->Jacobian(q);
    const FullMatrix omega_dot_gradient_squared = OmegaDotGradientSquared(
        q, quadrature::QuadraturePointIndex(angle_index));
    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) +=
            inverse_sigma_t * omega_dot_gradient_squared(i,j) * jacobian;
      }
    }
  }
}

template <int dim>
std::vector<double> CFEMSelfAdjointAngularFlux<dim>::OmegaDotGradient(
    int cell_quadrature_point,
    quadrature::QuadraturePointIndex angular_index) const {
  std::vector<double> return_vector(cell_degrees_of_freedom_);
  std::pair<int, int> index{cell_quadrature_point, angular_index.get()};
  for (int i = 0; i < cell_degrees_of_freedom_; ++i)
    return_vector.at(i) = omega_dot_gradient_.at(index)[i];
  return return_vector;
}

template <int dim>
FullMatrix CFEMSelfAdjointAngularFlux<dim>::OmegaDotGradientSquared(
    int cell_quadrature_point,
    quadrature::QuadraturePointIndex angular_index) const {
  return omega_dot_gradient_squared_.at(
      {cell_quadrature_point, angular_index.get()});
}

// PRIVATE FUNCTIONS ===========================================================
template <int dim>
void CFEMSelfAdjointAngularFlux<dim>::ValidateAndSetCell(
    const bart::formulation::CellPtr<dim> &cell_ptr,
    std::string function_name) {
  std::string error{"Error in CFEMSelfAdjointAngularFlux function " +
      function_name + ": passed cell pointer is invalid"};
  AssertThrow(cell_ptr.state() == dealii::IteratorState::valid,
              dealii::ExcMessage(error))
  finite_element_ptr_->SetCell(cell_ptr);
}

template <int dim>
void CFEMSelfAdjointAngularFlux<dim>::ValidateMatrixSize(
    const bart::formulation::FullMatrix& to_validate,
    std::string called_function_name) {
  auto [rows, cols] = std::pair{to_validate.n_rows(), to_validate.n_cols()};

  std::ostringstream error_string;
  error_string << "Error in CFEMSelfAdjointAngularFlux function "
               << called_function_name
               <<": passed matrix size is invalid, expected size ("
               << cell_degrees_of_freedom_ << ", " << cell_degrees_of_freedom_
               << "), actual size: (" << rows << ", " << cols << ")";

  AssertThrow((static_cast<int>(rows) == cell_degrees_of_freedom_) &&
      (static_cast<int>(cols) == cell_degrees_of_freedom_),
      dealii::ExcMessage(error_string.str()))
}

template <int dim>
void CFEMSelfAdjointAngularFlux<dim>::ValidateVectorSize(
    const bart::formulation::Vector& to_validate,
    std::string called_function_name) {

  int rows = to_validate.size();

  std::ostringstream error_string;
  error_string << "Error in CFEMSelfAdjointAngularFlux function "
               << called_function_name
               <<": passed vector size is invalid, expected size ("
               << cell_degrees_of_freedom_ << ", 1), actual size: (" << rows
               << ", 1)";

  AssertThrow((static_cast<int>(rows) == cell_degrees_of_freedom_),
              dealii::ExcMessage(error_string.str()))
}

template class CFEMSelfAdjointAngularFlux<1>;
template class CFEMSelfAdjointAngularFlux<2>;
template class CFEMSelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
