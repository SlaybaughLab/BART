#ifndef BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_

#include "calculator/residual/cell_isotropic_residual_i.hpp"

#include "domain/finite_element/finite_element_i.hpp"
#include "data/cross_sections/cross_sections_i.hpp"

#include "utility/has_dependencies.h"

namespace bart::calculator::residual {

template <int dim>
class CellIsotropicResidual : public CellIsotropicResidualI<dim>, public utility::HasDependencies {
 public:
  using typename CellIsotropicResidualI<dim>::CellPtr;
  using typename CellIsotropicResidualI<dim>::FluxMoments;
  using typename CellIsotropicResidualI<dim>::Vector;

  using CrossSections = data::cross_sections::CrossSectionsI;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;

  CellIsotropicResidual(std::shared_ptr<CrossSections>, std::shared_ptr<FiniteElement>);
  auto CalculateCellResidual(CellPtr, FluxMoments*, FluxMoments*, int group) -> Vector override;

  auto cross_sections_ptr() { return cross_sections_ptr_.get(); }
  auto finite_element_ptr() { return finite_element_ptr_.get(); }
 private:
  std::shared_ptr<CrossSections> cross_sections_ptr_{ nullptr };
  std::shared_ptr<FiniteElement> finite_element_ptr_{ nullptr };
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_
