#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_

#include "system/moments/spherical_harmonic_types.h"

namespace bart {

namespace system {

namespace solution {

class MPIAngularI;

} // namespace solution

} // namespace system

namespace quadrature {

namespace angular {

template <int dim>
class AngularQuadratureSetI;

} // namespace angular

namespace calculators {

template <int dim>
class SphericalHarmonicMomentsI {
 public:
  virtual ~SphericalHarmonicMomentsI() = default;

  virtual system::moments::MomentVector CalculateMoment(
      system::solution::MPIAngularI* solution,
      system::moments::HarmonicL harmonic_l,
      system::moments::HarmonicL harmonic_m) const = 0;

  virtual angular::AngularQuadratureSetI<dim>* angular_quadrature_set_ptr() const = 0;
};

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_