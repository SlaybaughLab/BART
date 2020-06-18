#ifndef BART_SRC_SYSTEM_SOLUTION_SOLUTION_TYPES_H_
#define BART_SRC_SYSTEM_SOLUTION_SOLUTION_TYPES_H_

#include <map>

#include "system/solution/mpi_group_angular_solution_i.h"
#include "system/system_types.h"

namespace bart {

namespace system {

namespace solution {

using AngularSolutionPtr = std::shared_ptr<dealii::Vector<double>>;

using EnergyGroupToAngularSolutionPtrMap = std::map<SolutionIndex, AngularSolutionPtr>;

} // namespace solution

} // namespace system

} // namespace bart

#endif //BART_SRC_SYSTEM_SOLUTION_SOLUTION_TYPES_H_
