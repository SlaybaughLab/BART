#ifndef BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_

#include "utility/has_description.h"

#include <deal.II/lac/vector.h>

#include "system/moments/spherical_harmonic_i.h"

namespace bart::calculator::residual {

/*! \brief Interface for classes that calculate the isotropic scattering residual for a whole domain.
 *
 * In general, these classes will use a CellIsotropicResidualI class to calculate the residual for each cell \f$K\f$ and
 * combine them into a Vector that gives the total isotropic scattering residual at each global degree of freedom summed
 * over all groups,
 * \f[
 * \big< R \big> = \sum_{K}^{T_K}\sum_{g = 0}^{G - 1}R_{K,g, 0}
 * \f]
 */
class DomainIsotropicResidualI : public utility::HasDescription {
 public:
  using Vector = dealii::Vector<double>;
  using FluxMoments = system::moments::SphericalHarmonicI;
  virtual ~DomainIsotropicResidualI() = default;

  /*! \brief Calculate the domain isotropic scattering residual.
   *
   * @param current_scalar_flux_ the scalar flux for step \f$k + 1/2\f$
   * @param previous_scalar_flux_ the scalar flux for step \f$k\f$
   */
  virtual auto CalculateDomainResidual(FluxMoments* current_flux_moments,
                                       FluxMoments* previous_flux_moments) -> Vector = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
