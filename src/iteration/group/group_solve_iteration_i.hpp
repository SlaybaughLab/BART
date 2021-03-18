#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_HPP_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_HPP_

#include "utility/has_description.h"
#include "system/solution/solution_types.h"
#include "system/system.hpp"

//! Classes that provide within-group iterations.
namespace bart::iteration::group {

/*! \brief Interface for classes that perform within-group iterations.
 *
 * The Iterate function should update a system to a point where all groups have converged. Most of the time, specific
 * angular solutions are not stored, so this interface also provides a method for giving a mapping of angular flux
 * solutions to update.
 *
 */
class GroupSolveIterationI : public utility::HasDescription {
 public:
  using EnergyGroupToAngularSolutionPtrMap = system::solution::EnergyGroupToAngularSolutionPtrMap;
  virtual ~GroupSolveIterationI() = default;
  /*! \brief Iterate system until all groups have converged. */
  virtual auto Iterate(system::System &system) -> void = 0;
  /*! \brief Provide a mapping of group to pointers to angular solutions to update. */
  virtual auto UpdateThisAngularSolutionMap(EnergyGroupToAngularSolutionPtrMap) -> GroupSolveIterationI& = 0;
};

} // namespace bart::iteration::group

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_I_HPP_
