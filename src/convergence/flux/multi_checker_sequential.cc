#include "multi_checker_sequential.h"

namespace bart {

namespace convergence {

namespace flux {

bool MultiCheckerSequential::CheckIfConverged(
    data::MultiFluxPtrs &current_iteration,
    data::MultiFluxPtrs &last_iteration) {
  
  AssertThrow(current.size() > 0,
              dealii::ExcMessage("Current iteration fluxes is empty"));
  AssertThrow(last.size() > 0,
              dealii::ExcMessage("Last iteration fluxes is empty"));
  AssertThrow(last.size() == current.size(),
              dealii::ExcMessage("Current & last group iteration fluxes have"
                                 "different sizes"));

  for (auto current_pair = current.cbegin(), last_pair = last.cbegin();
       current_pair != current.cend(), last_pair != last.cend();
       ++current_pair, ++last_pair) {

    auto &[current_index, current_group_flux_ptr] = *current_pair;
    auto &[last_index, last_group_flux_ptr] = *last_pair;
    
    AssertThrow(current_index == last_index,
                dealii::ExcMessage("Current & last group indices mismatched"));

    if (!tester_->CheckIfConverged(*current_group_flux_ptr,
                                   *last_group_flux_ptr)) {
      failed_index_ = current_index;
      return false;
    }
  }
  return true;
}

} // namespace flux

} // namespace convergence

} // namespace bart
