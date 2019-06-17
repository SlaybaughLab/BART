#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_

#include <memory>

#include "quadrature/angular/angular_quadrature_set_i.h"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template <int dim>
class SphericalHarmonicMoments : public SphericalHarmonicMomentsI {
 public:
  virtual ~SphericalHarmonicMoments() = default;

 protected:
  std::shared_ptr<angular::AngularQuadratureSetI<dim>> angular_quadrature_ptr_;
};


} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_