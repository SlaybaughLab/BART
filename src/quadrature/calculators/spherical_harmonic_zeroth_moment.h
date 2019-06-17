#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT

#include "quadrature/calculators/spherical_harmonic_moments.h"

namespace bart {

namespace quadrature {

namespace calculators {

template <int dim>
class SphericalHarmonicZerothMoment : public SphericalHarmonicMoments<dim> {
 public:
  SphericalHarmonicZerothMoment(
      std::shared_ptr<angular::AngularQuadratureSetI<dim>> angular_quadrature_ptr)
      : SphericalHarmonicMoments<dim>(angular_quadrature_ptr) {}

  virtual ~SphericalHarmonicZerothMoment() = default;
};

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT