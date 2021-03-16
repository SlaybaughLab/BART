#ifndef BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
#define BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_

#include "utility/has_description.h"

#include <deal.II/lac/vector.h>

#include "system/moments/spherical_harmonic_i.h"

namespace bart::calculator::residual {

class DomainIsotropicResidualI : public utility::HasDescription {
 public:
  using Vector = dealii::Vector<double>;
  using FluxMoments = system::moments::SphericalHarmonicI;

  virtual ~DomainIsotropicResidualI() = default;

  virtual auto CalculateDomainResidual(FluxMoments* current_flux_moments,
                                       FluxMoments* previous_flux_moments) -> Vector = 0;
};

} // namespace bart::calculator::residual

#endif //BART_SRC_CALCULATOR_RESIDUAL_DOMAIN_ISOTROPIC_RESIDUAL_I_HPP_
