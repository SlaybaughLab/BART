#ifndef BART_SRC_CONVERGENCE_FLUX_L1_THRESHOLD_H_
#define BART_SRC_CONVERGENCE_FLUX_L1_THRESHOLD_H_

#include "flux_i.h"

namespace bart {

namespace convergence {

/*! \brief Checks for convergence between two provided fluxes. */

class FluxL1Threshold : public FluxI {
  FluxL1Threshold() = default;
  ~FluxL1Threshold() = default;
  bool isConverged(data::Flux &, data::Flux &) override;
};
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_L1_THRESHOLD_H_
