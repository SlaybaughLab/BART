#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_

#include "system/solution/mpi_angular_i.h"

namespace bart {

namespace system {

namespace solution {

class MPIAngular : public MPIAngularI {
 public:

  MPIAngular(const int total_groups, const int total_angles = 1);

  int total_angles() const override { return total_angles_; }
  int total_groups() const override { return total_groups_; }
  const SolutionMap& solutions() const override { return solutions_; };

 private:
  SolutionMap solutions_;
  const int total_angles_;
  const int total_groups_;
};

} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_