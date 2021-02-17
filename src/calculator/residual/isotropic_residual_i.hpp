#ifndef BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_I_HPP_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include "system/moments/spherical_harmonic_i.h"

namespace bart::calculator::residual {

/*! \brief Interface for classes that calculate the isotropic component of the residual.
 *
 * The isotropic residual for group g is defined in Adams & Morel (1993) as
 * \f[
 * R^{k + 1/2}_{g, 0} = \sum_{g' = g + 1}^G \Sigma_{s, g' \to g, 0}\left(\phi_{g', 0}^{k + 1/2} - \phi_{g', 0}^{k}\right)
 * \f]
 *
*/
class IsotropicResidualI {
 public:
  using FullMatrix = dealii::FullMatrix<double>;
  using Moments = system::moments::SphericalHarmonicI;
  using Vector = dealii::Vector<double>;
  virtual ~IsotropicResidualI() = default;
  /*! \brief Calculate the isotropic residual
   *
   * @param half_step_scalar_flux_ptr the flux moments for step \f$k + 1/2\f$
   * @param previous_step_scalar_flux_ptr the flux moments for step \f$k\f$
   * @param group group to calculate, \f$g\f$
   * @param sigma_s scattering cross-section
   * @return a vector containing the calculated isotropic residual
   */
  virtual auto CalculateIsotropicResidual(Moments* half_step_scalar_flux_ptr, Moments* previous_step_scalar_flux_ptr,
                                          const int group, const FullMatrix& sigma_s) const -> Vector = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_ISOTROPIC_RESIDUAL_I_HPP_
