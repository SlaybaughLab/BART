#include "formulation/equation/diffusion.h"
#include "diffusion.h"

namespace bart {

namespace formulation {

namespace equation {

template<int dim>
void Diffusion<dim>::Precalculate(const CellPtr &cell_ptr) {
  SetCell(cell_ptr);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    Matrix shape_squared(cell_degrees_of_freedom_,
                         cell_degrees_of_freedom_);
    Matrix gradient_squared(cell_degrees_of_freedom_,
                            cell_degrees_of_freedom_);

    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i,j) = finite_element_->values()->shape_value(i,q) *
            finite_element_->values()->shape_value(j,q);
        gradient_squared(i,j) = finite_element_->values()->shape_grad(i,q) *
            finite_element_->values()->shape_grad(j,q);
      }
    }

    shape_squared_[q] = shape_squared;
    gradient_squared_[q] = gradient_squared;
  }
}

template<int dim>
void Diffusion<dim>::FillCellBilinearTerm(Matrix &to_fill,
                                          const CellPtr &cell_ptr,
                                          const GroupNumber group) const {
  SetCell(cell_ptr);
  MaterialID material_id = cell_ptr->material_id();

  double sigma_t = cross_sections_->sigma_t.at(material_id)[group];
  double sigma_s = cross_sections_->sigma_s.at(material_id)(group, group);
  double sigma_r = sigma_t - sigma_s;

  double diff_coef = cross_sections_->diffusion_coef.at(material_id)[group];

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    double gradient = finite_element_->values()->JxW(q);
    to_fill.add(diff_coef*gradient, gradient_squared_[q]);
    to_fill.add(sigma_r*gradient, shape_squared_[q]);
  }
}



template <int dim>
void Diffusion<dim>::FillCellLinearScatteringTerm(Matrix &to_fill,
                                                  const CellPtr &cell_ptr,
                                                  const GroupNumber group) const {
  SetCell(cell_ptr);
  MaterialID material_id = cell_ptr->material_id();



}
template<int dim>
void Diffusion<dim>::FillBoundaryBilinearTerm(
    Matrix &to_fill,
    const CellPtr &cell_ptr,
    const FaceNumber face_number,
    const BoundaryType boundary_type) const {

  if (boundary_type == BoundaryType::kVacuum) {
    SetFace(cell_ptr, face_number);

    for (int q = 0; q < face_quadrature_points_; ++q) {
      for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
        for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
          to_fill(i,j) += (
              finite_element_->face_values()->shape_value(i, q) *
              finite_element_->face_values()->shape_value(j,q)
              ) * 0.5 * finite_element_->face_values()->JxW(q);
        }
      }
    }
  }
}

template class Diffusion<1>;
template class Diffusion<2>;
template class Diffusion<3>;

} // namespace equation

} // namespace formulation

} // namespace bart