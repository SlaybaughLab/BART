#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "system/solution/mpi_group_angular_solution_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template<int dim>
system::moments::MomentVector SphericalHarmonicZerothMoment<dim>::CalculateMoment(
    system::solution::MPIGroupAngularSolutionI* solution,
    system::GroupNumber group,
    system::moments::HarmonicL,
    system::moments::HarmonicL) const {

  // Verify that the solution and the angular quadrature set have the same
  // number of angles.

  const int total_angles = solution->total_angles();

  AssertThrow(
      angular_quadrature_ptr_->total_quadrature_points() == total_angles,
      dealii::ExcMessage("Error: angular quadrature set and solution must "
                         "have the same number of angles."));

  auto weights = angular_quadrature_ptr_->quadrature_weights();

  system::moments::MomentVector return_vector;

  for (int angle = 0; angle < total_angles; ++angle) {
    auto mpi_solution = (*solution)[{group, angle}];
    system::moments::MomentVector angle_vector(mpi_solution);

    if (return_vector.size() == 0) {
      return_vector.reinit(angle_vector);
      return_vector = 0.0;
    }

    return_vector.add(weights[angle], angle_vector);
  }

  return return_vector;
}

template class SphericalHarmonicZerothMoment<1>;
template class SphericalHarmonicZerothMoment<2>;
template class SphericalHarmonicZerothMoment<3>;

} // namespace calculators

} // namespace quadrature

} // namespace bart