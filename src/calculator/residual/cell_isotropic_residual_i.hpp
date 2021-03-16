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
 * This class calculates the integrated isotropic scattering residual for a specified cell \f$K \in T_K\f$ based on
 * provided scalar fluxes and for a specified group using the cell quadrature:
 *
 * * \f[
 * R_{K,g, 0} = \sum_{q = 0}^Q\left[\sum_{g' = g + 1}^G \sigma_{s, g' \to g}\left(\phi^{k + 1/2}_{K,g', 0}(q) - \phi^{k}_{K,g', 0}(q)\right)\right]J_K(q)
 * \f]
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

  ~CellIsotropicResidualI() = default;
  /*! \brief Calculate the cell residual.
   *
   * Integrated and return the cell isotropic residual component as described.
   *
   * @param current_scalar_flux_ the scalar flux for step \f$k + 1/2\f$
   * @param previous_scalar_flux_ the scalar flux for step \f$k\f$
   * @param group the group to calculate the residual for
   * @return double value of the cell isotropic scattering residual
   */
  virtual auto CalculateCellResidual(CellPtr, FluxMoments* current_scalar_flux_, FluxMoments* previous_scalar_flux_,
                                     int group) -> double = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_CELL_ISOTROPIC_RESIDUAL_I_HPP_
