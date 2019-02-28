#include "convergence/final_flux.h"

namespace bart {

namespace convergence {

Status FinalFlux::CheckFinalConvergence() {
  return convergence_status_;
}

void FinalFlux::ProvideMultiChecker(
    std::unique_ptr<flux::MultiCheckerI> &checker) {
  checker_ = std::move(checker);
}

} // namespace convergence

} // namespace bart
