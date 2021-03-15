#include "calculator/residual/cell_isotropic_residual.hpp"

namespace bart::calculator::residual {

template<int dim>
CellIsotropicResidual<dim>::CellIsotropicResidual(std::shared_ptr<CrossSections> cross_sections_ptr,
                                                  std::shared_ptr<FiniteElement> finite_element_ptr)
    : cross_sections_ptr_(std::move(cross_sections_ptr)),
      finite_element_ptr_(std::move(finite_element_ptr)) {
  std::string function_name{ "CellIsotropicResidual constructor"};
  this->AssertPointerNotNull(cross_sections_ptr_.get(), "cross-sections", function_name);
  this->AssertPointerNotNull(finite_element_ptr_.get(), "finite element", function_name);
}

template<int dim>
auto CellIsotropicResidual<dim>::CalculateCellResidual(CellPtr, FluxMoments*, FluxMoments*, const int) -> Vector {
  return Vector();
}

template class CellIsotropicResidual<1>;
template class CellIsotropicResidual<2>;
template class CellIsotropicResidual<3>;

} // namespace bart::calculator::residual
