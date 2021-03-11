#include "formulation/scalar/two_grid_diffusion.hpp"

namespace bart::formulation::scalar {

template<int dim>
TwoGridDiffusion<dim>::TwoGridDiffusion(std::shared_ptr<FiniteElement> finite_element_ptr,
                                        std::shared_ptr<CrossSections> cross_sections_ptr,
                                        std::shared_ptr<OneGroupCrossSections> one_group_cross_sections_ptr)
    : one_group_cross_sections_ptr_(std::move(one_group_cross_sections_ptr)),
      Diffusion<dim>(std::move(finite_element_ptr), std::move(cross_sections_ptr)) {
  this->template AssertPointerNotNull(one_group_cross_sections_ptr_.get(), "one group cross sections",
                                      "two grid diffusion constructor");
}
template<int dim>
auto TwoGridDiffusion<dim>::FillCellCollisionTerm(Matrix& to_fill,
                                                  const CellPtr& cell_ptr,
                                                  GroupNumber) const -> void {
  this->VerifyInitialized(__FUNCTION__);
  this->finite_element_ptr_->SetCell(cell_ptr);

  const double sigma_removal{ one_group_cross_sections_ptr_->SigmaRemoval(cell_ptr->material_id()) };
  for (int q = 0; q < this->cell_quadrature_points_; ++q) {
    const double constant{ sigma_removal * this->finite_element_ptr_->Jacobian(q) };
    for (int i = 0; i < this->cell_degrees_of_freedom_; ++i) {
      for (int j = 0; j < this->cell_degrees_of_freedom_; ++j) {
        to_fill(i, j) += constant * this->shape_squared_.at(q)(i, j);
      }
    }
  }
}

template class TwoGridDiffusion<1>;
template class TwoGridDiffusion<2>;
template class TwoGridDiffusion<3>;

} // namespace bart::formulation::scalar
