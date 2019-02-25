#include "group_flux_checker_sequential.h"

namespace bart {

namespace convergence {

GroupFluxCheckerSequential::GroupFluxCheckerSequential(
    std::unique_ptr<FluxCheckerI> &tester)
    : tester_(std::move(tester)) {}

bool GroupFluxCheckerSequential::CheckIfConverged(data::ScalarGroupFluxPtrs &current,
                                                  data::ScalarGroupFluxPtrs &last) {
  AssertThrow(current.size() > 0,
              dealii::ExcMessage("Current iteration group fluxes has 0 groups"));
  AssertThrow(last.size() > 0,
              dealii::ExcMessage("Last iteration group fluxes has 0 groups"));
  AssertThrow(last.size() == current.size(),
              dealii::ExcMessage("Current & last group iteration fluxes have"
                                 "different sizes"));

  for (auto current_pair = current.cbegin(), last_pair = last.cbegin();
       current_pair != current.cend(), last_pair != last.cend();
       ++current_pair, ++last_pair) {

    auto &[current_group_number, current_group_flux_ptr] = *current_pair;
    auto &[last_group_number, last_group_flux_ptr] = *last_pair;
    
    AssertThrow(current_group_number == last_group_number,
                dealii::ExcMessage("Current & last group numbers mismatched"));

    if (!tester_->CheckIfConverged(*current_group_flux_ptr,
                                   *last_group_flux_ptr)) {
      failed_group_ = current_group_number;
      return false;
    }
  }
  return true;
}

bool GroupFluxCheckerSequential::CheckIfConverged(data::AngularGroupFluxPtrs &current,
                                                  data::AngularGroupFluxPtrs &last) {
  return true;
}

} // namespace convergence

} // namespace bart
