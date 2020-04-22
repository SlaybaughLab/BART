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
      std::shared_ptr<QuadratureSetI<dim>> quadrature_set_ptr)
      : SphericalHarmonicMoments<dim>(quadrature_set_ptr) {}

  system::moments::MomentVector CalculateMoment(
      system::solution::MPIGroupAngularSolutionI* solution,
      system::GroupNumber group,
      system::moments::HarmonicL harmonic_l,
      system::moments::HarmonicL harmonic_m) const override;

  virtual ~SphericalHarmonicZerothMoment() = default;

 protected:
  using SphericalHarmonicMoments<dim>::quadrature_set_ptr_;
};

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif // BART_SRC_QUADRATURE_CALCULATORS_SPHERICAL_HARMONIC_ZEROTH_MOMENT