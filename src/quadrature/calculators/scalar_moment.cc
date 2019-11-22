#include "quadrature/calculators/scalar_moment.h"

namespace bart {

namespace quadrature {

namespace calculators {

system::moments::MomentVector ScalarMoment::CalculateMoment(
    system::solution::MPIGroupAngularSolutionI *solution,
    system::GroupNumber /*group*/,
    system::moments::HarmonicL /*harmonic_l*/,
    system::moments::HarmonicL /*harmonic_m*/) const {

  system::moments::MomentVector return_vector;

  return return_vector;
}

} // namespace calculators

} // namespace quadrature

} // namespace bart

