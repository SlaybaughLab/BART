#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "system/solution/mpi_group_angular_solution_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

template<int dim>
system::moments::MomentVector SphericalHarmonicZerothMoment<dim>::CalculateMoment(
    system::solution::MPIGroupAngularSolutionI* solution,
    system::GroupNumber,
    system::moments::HarmonicL,
    system::moments::HarmonicL) const {

  // Verify that the solution and the angular quadrature set have the same
  // number of angles.

  const int total_angles = solution->total_angles();
  const int quadrature_size = this->quadrature_set_ptr_->size();

  AssertThrow(quadrature_size == total_angles,
      dealii::ExcMessage("Error: angular quadrature set and solution must "
                         "have the same number of angles."))

  system::moments::MomentVector return_vector;

  for (auto quadrature_point_ptr : *quadrature_set_ptr_) {
    const int angle_index =
        quadrature_set_ptr_->GetQuadraturePointIndex(quadrature_point_ptr);
    auto mpi_solution = solution->GetSolution(angle_index);

    system::moments::MomentVector angle_vector(mpi_solution);
    const double quadrature_point_weight = quadrature_point_ptr->weight();

    if (return_vector.size() == 0) {
      return_vector.reinit(angle_vector);
      return_vector = 0.0;
    }

    return_vector.add(quadrature_point_weight, angle_vector);

  }

  return return_vector;
}

template class SphericalHarmonicZerothMoment<1>;
template class SphericalHarmonicZerothMoment<2>;
template class SphericalHarmonicZerothMoment<3>;

} // namespace calculators

} // namespace quadrature

} // namespace bart