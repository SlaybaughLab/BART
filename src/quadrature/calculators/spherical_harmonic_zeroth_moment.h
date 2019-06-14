#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT

#include "quadrature/calculators/spherical_harmonic_moments_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

class SphericalHarmonicZerothMoment : public SphericalHarmonicMomentsI {
 public:
  virtual ~SphericalHarmonicZerothMoment() = default;
};

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT