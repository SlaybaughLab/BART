#include "group_flux_checker_sequential.h"

namespace bart {

namespace convergence {

GroupFluxCheckerSequential::GroupFluxCheckerSequential(
    std::unique_ptr<FluxCheckerI> &tester)
    : tester_(std::move(tester)) {}

bool GroupFluxCheckerSequential::isConverged(data::GroupFluxPointers &current,
                                             data::GroupFluxPointers &last) {
  AssertThrow(current.size() > 0,
              dealii::ExcMessage("Current iteration group fluxes has 0 groups"));
  AssertThrow(last.size() > 0,
              dealii::ExcMessage("Last iteration group fluxes has 0 groups"));
  AssertThrow(last.size() == current.size(),
              dealii::ExcMessage("Current & last group iteration fluxes have"
                                 "different sizes"));

  for (auto current_flux = current.cbegin(), last_flux = last.cbegin();
       current_flux != current.cend(), last_flux != last.cend();
       ++current_flux, ++last_flux) {
    
    AssertThrow(std::get<const data::Group>(*current_flux) == std::get<const data::Group>(*last_flux),
                dealii::ExcMessage("Current & last group numbers mismatched"));
    
    if (!tester_->isConverged(*(std::get<1>(*current_flux)),
                              *(std::get<1>(*last_flux)))) {
      failed_group_ = std::get<const data::Group>(*current_flux);
      converged_ = false;
      return false;
    }
  }
  converged_ = true;
  return true;
}

} // namespace convergence

} // namespace bart
