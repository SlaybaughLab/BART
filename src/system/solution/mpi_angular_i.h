#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_

#include <map>

#include "data/system/system_types.h"
#include "system/solution/solution_i.h"

namespace bart {

namespace system {

namespace solution {

class MPIAngularI :
    public SolutionI<data::system::Index, data::system::MPIVector> {
 public:
  using SolutionMap = std::map<data::system::Index, data::system::MPIVector>;
  using MPIVector = data::system::MPIVector;
  using Index = data::system::Index;

  virtual int total_groups() const = 0;
  virtual int total_angles() const = 0;
};


} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_