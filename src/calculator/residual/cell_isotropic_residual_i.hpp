#ifndef BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_

#include "domain/domain_types.hpp"
#include "system/moments/spherical_harmonic_i.h"

//! Calculators for determining residuals
namespace bart::calculator::residual {

/*! \brief Interface for classes that calculate the isotropic scattering residual as defined in the two-grid method.
 *
 * For the two-grid method of <a href="https://doi.org/10.13182/NSE115-253">Adams and Morel (1993)</a>, the isotropic
 * scattering component of the residual is,
 * \f[
 * R_{g, 0} = \sum_{g' = g + 1}^G \sigma_{s, g' \to g}\left(\phi^{k + 1/2}_{g', 0} - \phi^{k}_{g', 0}\right)
 * \f]
 *
 * This class calculates the isotropic scattering residual for a specified cell \f$K \in T_K\f$ based on
 * provided scalar fluxes and for a specified group.
 *
 *
 * @tparam dim spatial dimension.
 */
template <int dim>
class CellIsotropicResidualI {
 public:
  //! Cell pointer type
  using CellPtr = typename domain::CellPtr<dim>;
  //! Flux moments type
  using FluxMoments = system::moments::SphericalHarmonicI;

  virtual ~CellIsotropicResidualI() = default;
  /*! \brief Calculate the cell residual.
   */
  virtual auto CalculateCellResidual(dealii::Vector<double>& to_fill, CellPtr, FluxMoments* current_scalar_flux_,
                                     FluxMoments* previous_scalar_flux_, int group) -> void = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_
