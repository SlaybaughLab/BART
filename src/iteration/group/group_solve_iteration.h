#ifndef BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
#define BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_

#include "iteration/group/group_solve_iteration_i.h"

namespace bart {

namespace iteration {

namespace group {

class GroupSolveIteration : public GroupSolveIterationI {
 public:
  virtual ~GroupSolveIteration() = default;
};

} // namespace group

} // namespace iteration



} //namespace bart

#endif //BART_SRC_ITERATION_GROUP_GROUP_SOLVE_ITERATION_H_
