#include "convergence/final_checker_or_n.h"

#include "data/moment_types.h"
#include "convergence/moments/single_moment_checker_i.h"
#include "convergence/moments/multi_moment_checker_i.h"

namespace bart {

namespace convergence {

template <>
Status FinalCheckerOrN<data::MomentVector ,
                       moments::SingleMomentCheckerI>::CheckFinalConvergence(
    data::MomentVector & current_iteration,
    data::MomentVector & previous_iteration) {

  StatusDeltaAndIterate(current_iteration, previous_iteration);
  return convergence_status_;
}

template <>
Status FinalCheckerOrN<data::MomentsMap ,
                       moments::MultiMomentCheckerI>::CheckFinalConvergence(
    data::MomentsMap & current_iteration,
    data::MomentsMap & previous_iteration) {

  StatusDeltaAndIterate(current_iteration, previous_iteration);
  convergence_status_.failed_index = checker_ptr_->failed_index();
  return convergence_status_;
}

template<typename CompareType, typename CheckerType>
void FinalCheckerOrN<CompareType, CheckerType>::StatusDeltaAndIterate(
    CompareType &current_iteration,
    CompareType &previous_iteration) {

  convergence_status_.is_complete =
      checker_ptr_->CheckIfConverged(current_iteration, previous_iteration);

  convergence_status_.delta = checker_ptr_->delta();

  ++convergence_status_.iteration_number;

  if (convergence_status_.iteration_number ==
      convergence_status_.max_iterations) {
    convergence_status_.is_complete = true;
  }
}

template class FinalCheckerOrN<data::MomentVector,
                               moments::SingleMomentCheckerI>;
template class FinalCheckerOrN<data::MomentsMap,
                               moments::MultiMomentCheckerI>;


} // namespace convergence

} // namespace bart