#include "convergence/final_checker_or_n.h"

#include "convergence/moments/single_moment_checker_i.h"

namespace bart {

namespace convergence {

template <typename CheckerType>
Status FinalCheckerOrN<CheckerType>::CheckFinalConvergence() {
  return convergence_status_;
}

template class FinalCheckerOrN<convergence::moments::SingleMomentCheckerI>;

} // namespace convergence

} // namespace bart