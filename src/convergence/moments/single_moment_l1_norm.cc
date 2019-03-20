#include "convergence/moments/single_moment_l1_norm.h"

namespace bart {

namespace convergence {

namespace moments {

bool SingleMomentL1Norm::CheckIfConverged(
    const data::MomentVector &current_iteration,
    const data::MomentVector &previous_iteration) {
//  data::FluxVector difference{current_iteration};
//  difference -= previous_iteration;
//  delta_ = difference.l1_norm()/current_iteration.l1_norm();
//  is_converged_ = delta_ <= max_delta_;
  return is_converged_;
}

} // namespace moments

} // namespace convergence

} // namespace bart
