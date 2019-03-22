#include "convergence/final_or_n.h"

#include "convergence/moments/single_moment_checker_i.h"

namespace bart {

namespace convergence {

template <typename CheckerType>
Status FinalOrN<CheckerType>::CheckFinalConvergence() {
  return convergence_status_;
}

template class FinalOrN<convergence::moments::SingleMomentCheckerI>;

} // namespace convergence

} // namespace bart