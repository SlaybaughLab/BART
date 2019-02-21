#include "group_flux_checker_sequential.h"

namespace bart {

namespace convergence {

GroupFluxCheckerSequential::GroupFluxCheckerSequential(
    std::unique_ptr<FluxCheckerI> &tester)
    : tester_(std::move(tester)) {}

bool GroupFluxCheckerSequential::isConverged(data::GroupFluxes &current,
                                             data::GroupFluxes &last) {
  AssertThrow(current.size() > 0,
              dealii::ExcMessage("Current iteration group fluxes has 0 groups"));
  AssertThrow(last.size() > 0,
              dealii::ExcMessage("Last iteration group fluxes has 0 groups"));
  AssertThrow(last.size() == current.size(),
              dealii::ExcMessage("Current & last group iteration fluxes have"
                                 "different sizes"));


};

} // namespace convergence

} // namespace bart
