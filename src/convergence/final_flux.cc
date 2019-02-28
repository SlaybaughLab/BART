#include "convergence/final_flux.h"

namespace bart {

namespace convergence {

Status FinalFlux::CheckFinalConvergence() {
  return convergence_status_;
}

} // namespace convergence

} // namespace bart
