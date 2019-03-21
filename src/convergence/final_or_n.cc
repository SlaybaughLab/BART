#include "convergence/final_or_n.h"

namespace bart {

namespace convergence {

Status FinalOrN::CheckFinalConvergence() {
  return convergence_status_;
}

} // namespace convergence

} // namespace bart