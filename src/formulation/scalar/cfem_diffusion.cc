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

template class CFEM_Diffusion<2>;

} // namespace scalar

} // namespace formulation

} // namespace bart