#ifndef BART_SRC_QUADRATURE_CALCULATORS_SCALAR_MOMENT_H_
#define BART_SRC_QUADRATURE_CALCULATORS_SCALAR_MOMENT_H_

#include "quadrature/calculators/spherical_harmonic_moments_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

class ScalarMoment : public SphericalHarmonicMomentsI {
 public:
  system::moments::MomentVector CalculateMoment(
      system::solution::MPIGroupAngularSolutionI *solution,
      system::GroupNumber group,
      system::moments::HarmonicL harmonic_l,
      system::moments::HarmonicL harmonic_m) const override;
};

} // namespace calculators

} // namespace quadrature

} //namespace bart

#endif //BART_SRC_QUADRATURE_CALCULATORS_SCALAR_MOMENT_H_
