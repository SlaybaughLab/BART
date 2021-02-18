#include "formulation/common/rhs_constant.hpp"

namespace bart::formulation::common {

template<int dim>
auto RHSConstant<dim>::FillCellConstantTerm(Vector& to_fill, const CellPtr& cell_ptr,
                                            const Vector& constant_vector) const -> void {
  finite_element_ptr_->SetCell(cell_ptr);
  const auto cell_degrees_of_freedom{ this->finite_element_ptr_->dofs_per_cell() };
  const auto cell_quadrature_points{ this->finite_element_ptr_->n_cell_quad_pts() };
  const auto constant_vector_at_quadrature{ this->finite_element_ptr_->ValueAtQuadrature(constant_vector) };

  for (int q = 0; q < cell_quadrature_points; ++q) {
    const double jacobian{ finite_element_ptr_->Jacobian(q) };
    for (int i = 0; i < cell_degrees_of_freedom; ++i) {
      to_fill(i) += constant_vector_at_quadrature.at(q) * finite_element_ptr_->ShapeValue(i, q) * jacobian;
    }
  }
}

template class RHSConstant<1>;
template class RHSConstant<2>;
template class RHSConstant<3>;

} // namespace bart::formulation::common
