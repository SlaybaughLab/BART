#include "convergence/final_flux_or_n.h"

namespace bart {

namespace convergence {

Status FinalFluxOrN::CheckFinalConvergence() {
  return convergence_status_;
}

} // namespace convergence

} // namespace bart
