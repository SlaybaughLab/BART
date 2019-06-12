#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_

#include "system/solution/mpi_angular_i.h"

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
class MPIAngular : public MPIAngularI {
 public:

  MPIAngular(const int total_groups, const int total_angles = 1);
  virtual ~MPIAngular() = default;

  int total_angles() const override { return total_angles_; }
  int total_groups() const override { return total_groups_; }
  const SolutionMap& solutions() const override { return solutions_; };

  /*! \brief Returns the system solution identified by the provided index */
  const MPIVector& operator[](const Index index) const override;
  /*! \brief Returns the system solution identified by the provided index */
        MPIVector& operator[](const Index index) override;

 private:
  SolutionMap solutions_;
  const int total_angles_;
  const int total_groups_;
};

} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_H_