#include "formulation/scalar/cfem_diffusion.h"

namespace bart {

namespace formulation {

namespace scalar {

template<int dim>
CFEM_Diffusion<dim>::CFEM_Diffusion(std::shared_ptr<domain::FiniteElementI<dim>> finite_element,
                                    std::shared_ptr<data::CrossSections> cross_sections)
    : finite_element_(finite_element),
      cross_sections_(cross_sections),
      cell_degrees_of_freedom_(finite_element->dofs_per_cell()),
      cell_quadrature_points_(finite_element->n_cell_quad_pts()),
      face_quadrature_points_(finite_element->n_face_quad_pts()) {}

template <int dim>
typename CFEM_Diffusion<dim>::InitializationToken
CFEM_Diffusion<dim>::Precalculate(const CellPtr cell_ptr) {

  finite_element_->SetCell(cell_ptr);

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    Matrix gradient_squared(cell_degrees_of_freedom_,
                            cell_degrees_of_freedom_);
    Matrix shape_squared(cell_degrees_of_freedom_,
                         cell_degrees_of_freedom_);

    for (int i = 0; i < cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < cell_degrees_of_freedom_; ++j) {
        shape_squared(i, j) =
            finite_element_->ShapeValue(i, q) *
            finite_element_->ShapeValue(j, q);
        gradient_squared(i, j) =
            finite_element_->ShapeGradient(i, q) *
            finite_element_->ShapeGradient(j, q);
      }
    }
    shape_squared_.push_back(shape_squared);
    gradient_squared_.push_back(gradient_squared);
  }
  InitializationToken token;
  return token;
}

template <int dim>
void CFEM_Diffusion<dim>::FillCellStreamingTerm(Matrix& to_fill,
                                                const InitializationToken,
                                                const CellPtr cell_ptr,
                                                const MaterialID material_id,
                                                const GroupNumber group) const {

  finite_element_->SetCell(cell_ptr);

  const double diffusion_coef =
      cross_sections_->diffusion_coef.at(material_id)[group];

  for (int q = 0; q < cell_quadrature_points_; ++q) {
    double jacobian = finite_element_->Jacobian(q);
    to_fill.add(diffusion_coef*jacobian, gradient_squared_[q]);
  }

}

template <int dim>
void CFEM_Diffusion<dim>::FillCellCollisionTerm(Matrix& to_fill,
                                                const InitializationToken,
                                                const CellPtr cell_ptr,
                                                const MaterialID material_id,
                                                const GroupNumber group) const {
}

template class CFEM_Diffusion<1>;
template class CFEM_Diffusion<2>;
template class CFEM_Diffusion<3>;

} // namespace scalar

} // namespace formulation

} // namespace bart