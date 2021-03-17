#ifndef BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_
#define BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_

#include <memory>

#include "quadrature/quadrature_set_i.hpp"
#include "quadrature/calculators/spherical_harmonic_moments_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template <int dim>
class SphericalHarmonicMoments : public SphericalHarmonicMomentsI {
 public:
  SphericalHarmonicMoments(
      std::shared_ptr<QuadratureSetI<dim>> quadrature_set_ptr)
  : quadrature_set_ptr_(quadrature_set_ptr) {}

  virtual ~SphericalHarmonicMoments() = default;

  QuadratureSetI<dim>* quadrature_set_ptr() const {
    return quadrature_set_ptr_.get();
  }

 protected:
  std::shared_ptr<QuadratureSetI<dim>> quadrature_set_ptr_;
};


} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_MOMENTS_H_