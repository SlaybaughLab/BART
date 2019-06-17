#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_

#include "quadrature/angular/angular_quadrature_set_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template <int dim>
class SphericalHarmonicMomentsI {
 public:
  virtual ~SphericalHarmonicMomentsI() = default;
  virtual angular::AngularQuadratureSetI<dim>* angular_quadrature_set_ptr() const = 0;
};

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_I_H_