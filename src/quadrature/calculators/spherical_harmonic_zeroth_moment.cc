#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "system/solution/mpi_angular_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template<int dim>
system::moments::MomentVector SphericalHarmonicZerothMoment<dim>::CalculateMoment(
    system::solution::MPIAngularI* solution,
    data::system::GroupNumber group,
    system::moments::HarmonicL,
    system::moments::HarmonicL) const {

  // Verify that the solution and the angular quadrature set have the same
  // number of angles.

  AssertThrow(
      angular_quadrature_ptr_->total_quadrature_points() == solution->total_angles(),
      dealii::ExcMessage("Error: angular quadrature set and solution must "
                         "have the same number of angles."));

  return bart::system::moments::MomentVector();
}

template class SphericalHarmonicZerothMoment<1>;
template class SphericalHarmonicZerothMoment<2>;
template class SphericalHarmonicZerothMoment<3>;

} // namespace calculators

} // namespace quadrature

} // namespace bart