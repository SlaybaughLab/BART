#include "convergence/moments/single_moment_checker_l1_norm.h"

namespace bart {

namespace convergence {

namespace moments {

bool SingleMomentCheckerL1Norm::IsConverged(
    const system::moments::MomentVector &current_iteration,
    const system::moments::MomentVector &previous_iteration) {
  system::moments::MomentVector difference(current_iteration);
  difference -= previous_iteration;
  delta_ = difference.l1_norm()/current_iteration.l1_norm();
  is_converged_ = delta_ <= max_delta_;
  return is_converged_;
}

} // namespace moments

} // namespace convergence

} // namespace bart
