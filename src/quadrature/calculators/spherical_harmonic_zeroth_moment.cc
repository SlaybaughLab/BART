#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"

namespace bart {

namespace quadrature {

namespace calculators {

template<int dim>
system::moments::MomentVector SphericalHarmonicZerothMoment<dim>::CalculateMoment(
    system::solution::MPIAngularI* solution,
    data::system::GroupNumber group,
    system::moments::HarmonicL harmonic_l,
    system::moments::HarmonicL harmonic_m) const {

  return bart::system::moments::MomentVector();
}

template class SphericalHarmonicZerothMoment<1>;
template class SphericalHarmonicZerothMoment<2>;
template class SphericalHarmonicZerothMoment<3>;

} // namespace calculators

} // namespace quadrature

} // namespace bart