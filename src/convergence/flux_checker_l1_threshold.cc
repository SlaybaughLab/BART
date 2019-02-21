#include "flux_checker_l1_threshold.h"

namespace bart {

namespace convergence {

bool FluxCheckerL1Threshold::isConverged(data::Flux &current_iteration,
                                         data::Flux &last_iteration) {
  data::Flux difference{current_iteration};
  difference -= last_iteration;
  return (difference.l1_norm()/current_iteration.l1_norm() <= threshold_);
}

} // namespace convergence

} // namespace bart
