#include "single_checker_delta_l1_norm.h"

namespace bart {

namespace convergence {

bool SingleCheckerDeltaL1Norm::CheckIfConverged(data::Flux &current_iteration,
                                                data::Flux &previous_iteration) {
  data::Flux difference{current_iteration};
  difference -= previous_iteration;
  delta_ = difference.l1_norm()/current_iteration.l1_norm();
  is_converged_ = delta_ <= max_delta_;
  return is_converged_;
}

} // namespace convergence

} // namespace bart
