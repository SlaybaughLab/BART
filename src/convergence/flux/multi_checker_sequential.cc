#include "multi_checker_sequential.h"

namespace bart {

namespace convergence {

namespace flux {

MultiCheckerSequential::MultiCheckerSequential(
    std::unique_ptr<SingleCheckerI> &checker) {
  ProvideChecker(checker);
}

bool MultiCheckerSequential::CheckIfConverged(
    data::MultiFluxPtrs &current_iteration,
    data::MultiFluxPtrs &previous_iteration) {
  
  AssertThrow(current_iteration.size() > 0,
              dealii::ExcMessage("Current iteration fluxes is empty"));
  AssertThrow(previous_iteration.size() > 0,
              dealii::ExcMessage("Previous iteration fluxes is empty"));
  AssertThrow(previous_iteration.size() == current_iteration.size(),
              dealii::ExcMessage("Current & previous group iteration fluxes have"
                                 "different sizes"));

  for (auto current_pair = current_iteration.cbegin(),
            previous_pair = previous_iteration.cbegin();
            current_pair != current_iteration.cend(),
            previous_pair != previous_iteration.cend();
            ++current_pair, ++previous_pair) {

    auto &[current_index, current_group_flux_ptr] = *current_pair;
    auto &[previous_index, previous_group_flux_ptr] = *previous_pair;
    
    AssertThrow(current_index == previous_index,
                dealii::ExcMessage("Current & previous group indices mismatched"));

    if (!checker_->CheckIfConverged(*current_group_flux_ptr,
                                    *previous_group_flux_ptr)) {
      failed_index_ = current_index;
      is_converged_ = false;
      failed_delta_ = checker_->delta();
      return false;
    }
  }
  failed_index_ = std::nullopt;
  failed_delta_ = std::nullopt;
  is_converged_ = true;
  return true;
}

} // namespace flux

} // namespace convergence

} // namespace bart
