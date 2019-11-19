#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_H_

#include "system/solution/mpi_group_angular_solution_i.h"

namespace bart {

namespace system {

namespace solution {

/*! \brief Creates and stores angular system solutions as PETSc MPI vectors.
 *
 * This default implementation instantiates the MPI Vectors during construction
 * based on the provided total number of angles and groups, which are assumed
 * to be constant.
 *
 */
class MPIGroupAngularSolution : public MPIGroupAngularSolutionI {
 public:

  MPIGroupAngularSolution(const int total_angles);
  virtual ~MPIGroupAngularSolution() = default;

  int total_angles() const override { return total_angles_; }

 private:
  const int total_angles_;
};

} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_H_