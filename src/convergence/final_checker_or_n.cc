#include "convergence/final_checker_or_n.h"

#include "data/moment_types.h"
#include "convergence/moments/single_moment_checker_i.h"

namespace bart {

namespace convergence {

template <typename CompareType, typename CheckerType>
Status FinalCheckerOrN<CompareType, CheckerType>::CheckFinalConvergence(
    CompareType& current_iteration,
    CompareType& previous_iteration
    ) {
  return convergence_status_;
}

template class FinalCheckerOrN<data::MomentVector,
                               convergence::moments::SingleMomentCheckerI>;


} // namespace convergence

} // namespace bart