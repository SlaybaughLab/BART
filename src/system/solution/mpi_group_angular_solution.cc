#include "system/solution/mpi_group_angular_solution.h"

namespace bart {

namespace system {

namespace solution {

MPIGroupAngularSolution::MPIGroupAngularSolution(const int total_angles)
    : total_angles_(total_angles) {
  for (int angle = 0; angle < total_angles; ++angle) {
    // Insert uninitilized MPIVectors
    MPIVector new_vector;
    solutions_[angle] = new_vector;
  }
  this->set_description("Group solution (MPI), total angles = " + total_angles);
  this->is_default_implementation(true);
}

} // namespace solution

} // namespace system

} // namespace bart