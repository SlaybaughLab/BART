#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_I_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_I_H_

#include <map>

#include "system/system_types.h"
#include "system/solution/solution.h"
#include "utility/has_description.h"

namespace bart {

namespace system {

namespace solution {

/*! \brief Interface for classes that store group angular system solutions as PETSc MPI vectors.
 *
 */
class MPIGroupAngularSolutionI :
    public Solution<system::AngleIndex, system::MPIVector> ,
    public utility::HasDescription {
 public:
  using SolutionMap = std::map<system::AngleIndex, system::MPIVector>;
  using MPIVector = system::MPIVector;
  using AngleIndex = system::AngleIndex;

  virtual ~MPIGroupAngularSolutionI() = default;

  /*! \brief Returns the total number of angles */
  virtual int total_angles() const = 0;
};


} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_MPI_GROUP_ANGULAR_SOLUTION_I_H_