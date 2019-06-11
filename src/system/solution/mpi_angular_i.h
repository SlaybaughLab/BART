#ifndef BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_
#define BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_

#include "data/system/system_types.h"
#include "system/solution/solution_i.h"

namespace bart {

namespace system {

namespace solution {

class MPIAngularI :
    public SolutionI<data::system::Index, data::system::MPIVector> {

};


} // namespace solution

} // namespace system

} // namespace bart

#endif // BART_SRC_SYSTEM_SOLUTION_MPI_ANGULAR_I_H_