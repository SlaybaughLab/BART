#ifndef BART_SRC_SOLVER_GMRES_H_
#define BART_SRC_SOLVER_GMRES_H_

#include "solver/linear_i.h"

namespace bart {

namespace solver {

class GMRES : public LinearI {
 public:
  GMRES() = default;
  ~GMRES() = default;
};

} // namespace solver

} // namespace bart

#endif // BART_SRC_SOLVER_GMRES_H_