#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_H_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_H_

#include "utility/has_description.h"
#include "system/solution/solution_types.h"

namespace bart {

namespace system {
class System;
} // namespace system

namespace iteration {

namespace group {

class GroupSolveIterationI : public utility::HasDescription {
 public:
  using EnergyGroupToAngularSolutionPtrMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  virtual ~GroupSolveIterationI() = default;
  virtual void Iterate(system::System &system) = 0;
  virtual auto UpdateThisAngularSolutionMap(EnergyGroupToAngularSolutionPtrMap&) -> GroupSolveIterationI& = 0;
};

} // namespace group

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_H_
