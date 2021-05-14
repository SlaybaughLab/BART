#ifndef BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_

#include "calculator/residual/cell_isotropic_residual_i.hpp"

#include "domain/finite_element/finite_element_i.hpp"
#include "data/cross_sections/cross_sections_i.hpp"

#include "utility/has_dependencies.h"

namespace bart::calculator::residual {

/*! \brief Default implementation of the cell isotropic scattering residual calculator. */
template <int dim>
class CellIsotropicResidual : public CellIsotropicResidualI<dim>, public utility::HasDependencies {
 public:
  using typename CellIsotropicResidualI<dim>::CellPtr;
  using typename CellIsotropicResidualI<dim>::FluxMoments;

  using CrossSections = data::cross_sections::CrossSectionsI;
  using FiniteElement = typename domain::finite_element::FiniteElementI<dim>;

  /*! \brief Constructor.
   * Takes shared ownership of cross-sections and finite-element objects.
   */
  CellIsotropicResidual(std::shared_ptr<CrossSections>, std::shared_ptr<FiniteElement>);

  auto CalculateCellResidual(dealii::Vector<double> &to_fill, CellPtr, FluxMoments* current_flux_moments_ptr,
                             FluxMoments* previous_flux_moments_ptr, int group) -> void override;

  /*! \brief Access cross-sections dependency. */
  auto cross_sections_ptr() { return cross_sections_ptr_.get(); }
  /*! \brief Access finite-element dependency. */
  auto finite_element_ptr() { return finite_element_ptr_.get(); }
 private:
  std::shared_ptr<CrossSections> cross_sections_ptr_{ nullptr }; //!< Cross-sections dependency for accessing \f$\sigma_s\f$
  std::shared_ptr<FiniteElement> finite_element_ptr_{ nullptr }; //!< Finite-element dependency for integration.
  const int n_cell_quadrature_points_; //!< total quadrature points per cell
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_HPP_
