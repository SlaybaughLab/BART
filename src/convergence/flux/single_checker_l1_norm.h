#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_

#include "single_checker_delta.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks for convergence between two fluxes using the percentage
 * change in the L1 norms */

class SingleCheckerDeltaL1Norm : public SingleCheckerDelta {
 public:
  SingleCheckerDeltaL1Norm() = default;
  ~SingleCheckerDeltaL1Norm() = default;
  bool CheckIfConverged(data::Flux &current_iteration,
                        data::Flux &previous_iteration) override;
};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_
