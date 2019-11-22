#include "quadrature/calculators/scalar_moment.h"

#include "system/solution/mpi_group_angular_solution_i.h"

namespace bart {

namespace quadrature {

namespace calculators {

system::moments::MomentVector ScalarMoment::CalculateMoment(
    system::solution::MPIGroupAngularSolutionI *solution,
    system::GroupNumber /*group*/,
    system::moments::HarmonicL /*harmonic_l*/,
    system::moments::HarmonicL /*harmonic_m*/) const {

  AssertThrow(solution->total_angles() == 1,
      dealii::ExcMessage("Error: Using ScalarMoment quadrature calculator "
                         "but solution appears to have more than one angle"));

  auto mpi_solution = solution->GetSolution(0);

  system::moments::MomentVector return_vector(mpi_solution);

  return return_vector;
}

} // namespace calculators

} // namespace quadrature

} // namespace bart

