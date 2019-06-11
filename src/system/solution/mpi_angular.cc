#include "system/solution/mpi_angular.h"

namespace bart {

namespace system {

namespace solution {

MPIAngular::MPIAngular(const int total_groups, const int total_angles)
    : total_angles_(total_angles),
      total_groups_(total_groups) {
  for (int group = 0; group < total_groups; ++group) {
    for (int angle = 0; angle < total_angles; ++angle) {
      // Insert uninitilized MPIVectors
      solutions_.insert(std::make_pair<Index, MPIVector>({group, angle}, {}));
    }
  }
};

} // namespace solution

} // namespace system

} // namespace bart