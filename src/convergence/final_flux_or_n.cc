#include "convergence/final_flux_or_n.h"

namespace bart {

namespace convergence {

Status FinalFluxOrN::CheckFinalConvergence() {
  return convergence_status_;
}
FinalFluxOrN::FinalFluxOrN(std::unique_ptr<flux::MultiCheckerI> checker,
                           std::shared_ptr<data::SystemFluxes> fluxes)
    : checker_(std::move(checker)),
      fluxes_(fluxes) {}

} // namespace convergence

} // namespace bart