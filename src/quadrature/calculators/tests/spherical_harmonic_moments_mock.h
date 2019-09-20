#ifndef BART_SRC_QUADRATURE_CALCULATORS_TESTS_SPHERICAL_HARMONIC_MOMENTS_MOCK_H_
#define BART_SRC_QUADRATURE_CALCULATORS_TESTS_SPHERICAL_HARMONIC_MOMENTS_MOCK_H_

#include "quadrature/calculators/spherical_harmonic_moments_i.h"

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace quadrature {

namespace calculators {

template <int dim>
class SphericalHarmonicMomentsMock : public SphericalHarmonicMomentsI<dim> {
 public:
  MOCK_CONST_METHOD4_T(CalculateMoment, system::moments::MomentVector(
          system::solution::MPIGroupAngularSolutionI* solution,
          system::GroupNumber group,
          system::moments::HarmonicL harmonic_l,
          system::moments::HarmonicL harmonic_m));
  MOCK_CONST_METHOD0_T(angular_quadrature_set_ptr,
      angular::AngularQuadratureSetI<dim>*());
  };

} // namespace calculators

} // namespace quadrature

} // namespace bart

#endif //BART_SRC_QUADRATURE_CALCULATORS_TESTS_SPHERICAL_HARMONIC_MOMENTS_MOCK_H_
