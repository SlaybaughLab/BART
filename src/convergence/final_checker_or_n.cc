#include "convergence/final_checker_or_n.h"

#include "system/moments/spherical_harmonic_types.h"
#include "convergence/moments/single_moment_checker_i.h"
#include "convergence/moments/multi_moment_checker_i.h"
#include "convergence/parameters/single_parameter_checker.hpp"

namespace bart {

namespace convergence {

template <typename CompareType, typename CheckerType>
Status FinalCheckerOrN<CompareType, CheckerType>::ConvergenceStatus(
    CompareType& current_iteration,
    CompareType& previous_iteration) {

  StatusDeltaAndIterate(current_iteration, previous_iteration);
  return convergence_status_;
}

template <>
Status FinalCheckerOrN<system::moments::MomentsMap ,
                       moments::MultiMomentCheckerI>::ConvergenceStatus(
    system::moments::MomentsMap & current_iteration,
    system::moments::MomentsMap & previous_iteration) {

  StatusDeltaAndIterate(current_iteration, previous_iteration);
  convergence_status_.failed_index = checker_ptr_->failed_index();
  return convergence_status_;
}

template<typename CompareType, typename CheckerType>
void FinalCheckerOrN<CompareType, CheckerType>::StatusDeltaAndIterate(
    CompareType &current_iteration,
    CompareType &previous_iteration) {

  convergence_status_.is_complete =
      checker_ptr_->IsConverged(current_iteration, previous_iteration);

  convergence_status_.delta = checker_ptr_->delta();

  ++convergence_status_.iteration_number;

  if (convergence_status_.iteration_number ==
      convergence_status_.max_iterations) {
    convergence_status_.is_complete = true;
  }
}

template class FinalCheckerOrN<system::moments::MomentVector,
                               moments::SingleMomentCheckerI>;
template class FinalCheckerOrN<system::moments::MomentsMap,
                               moments::MultiMomentCheckerI>;
template class FinalCheckerOrN<const system::moments::MomentsMap,
                               moments::MultiMomentCheckerI>;
template class FinalCheckerOrN<double, parameters::SingleParameterChecker>;


} // namespace convergence

} // namespace bart