#include "convergence/final_flux_or_n.h"

namespace bart {

namespace convergence {

FinalFluxOrN::FinalFluxOrN(std::unique_ptr<flux::MultiCheckerI> checker,
                           std::shared_ptr<data::SystemFluxes> fluxes)
    : checker_(std::move(checker)),
      fluxes_(fluxes) {}

Status FinalFluxOrN::CheckFinalConvergence() {
  convergence_status_.iteration_number += 1;
  bool flux_converged = checker_->CheckIfConverged(fluxes_->current_iteration,
                                                   fluxes_->previous_iteration);
  bool max_iterations_reached = (convergence_status_.iteration_number ==
      convergence_status_.max_iterations);

  convergence_status_.is_complete = (flux_converged || max_iterations_reached);
  convergence_status_.delta = checker_->delta();
  convergence_status_.failed_index = checker_->failed_index();

  return convergence_status_;
}

} // namespace convergence

} // namespace bart