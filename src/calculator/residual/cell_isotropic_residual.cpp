#include "calculator/residual/cell_isotropic_residual.hpp"

namespace bart::calculator::residual {

namespace  {
template <int dim>
auto GetCellQuadraturePoints(domain::finite_element::FiniteElementI<dim>* finite_element_ptr) -> int {
  AssertThrow(finite_element_ptr != nullptr,
              dealii::ExcMessage("Error in constructor in IntegratedFissionSource: finite_element_ptr is null"))
  return finite_element_ptr->n_cell_quad_pts();
}
} // namespace

template<int dim>
CellIsotropicResidual<dim>::CellIsotropicResidual(std::shared_ptr<CrossSections> cross_sections_ptr,
                                                  std::shared_ptr<FiniteElement> finite_element_ptr)
    : cross_sections_ptr_(std::move(cross_sections_ptr)),
      finite_element_ptr_(std::move(finite_element_ptr)),
      n_cell_quadrature_points_(GetCellQuadraturePoints(finite_element_ptr_.get())) {
  std::string function_name{ "CellIsotropicResidual constructor"};
  this->AssertPointerNotNull(cross_sections_ptr_.get(), "cross-sections", function_name);
  //this->AssertPointerNotNull(finite_element_ptr_.get(), "finite element", function_name);
}

template<int dim>
auto CellIsotropicResidual<dim>::CalculateCellResidual(dealii::Vector<double> &to_fill,
                                                       CellPtr cell_ptr,
                                                       FluxMoments* current_flux_moments_ptr,
                                                       FluxMoments* previous_flux_moments_ptr,
                                                       const int group) -> void {
  //finite_element_ptr_->SetCell(cell_ptr);
  const int total_groups = current_flux_moments_ptr->total_groups();
  const auto sigma_s{ cross_sections_ptr_->sigma_s().at(cell_ptr->material_id()) };
  const int cell_dofs = this->finite_element_ptr_->dofs_per_cell();
  std::vector<unsigned int> cell_global_dofs_indices(cell_dofs);
  cell_ptr->get_dof_indices(cell_global_dofs_indices);

  for (int group_in = group + 1; group_in < total_groups; ++group_in) {
    auto current_flux = current_flux_moments_ptr->GetMoment({group_in, 0, 0});
    auto previous_flux = previous_flux_moments_ptr->GetMoment({group_in, 0, 0});
    for (int i = 0; i < cell_dofs; ++i) {
      const auto index = cell_global_dofs_indices.at(i);
      to_fill(index) += sigma_s(group, group_in) * (current_flux(index) - previous_flux(index));
    }
  }
}

template class CellIsotropicResidual<1>;
template class CellIsotropicResidual<2>;
template class CellIsotropicResidual<3>;

} // namespace bart::calculator::residual
