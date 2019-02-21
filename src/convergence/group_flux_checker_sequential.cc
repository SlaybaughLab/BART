#include "group_flux_checker_sequential.h"

namespace bart {

namespace convergence {

GroupFluxCheckerSequential::GroupFluxCheckerSequential(
    std::unique_ptr<FluxCheckerI> &tester)
    : tester_(std::move(tester)) {}

} // namespace convergence

} // namespace bart
