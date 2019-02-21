#include "group_fluxes_sequential.h"

namespace bart {

namespace convergence {

GroupFluxesSequential::GroupFluxesSequential(std::unique_ptr<FluxI> &tester)
    : tester_(std::move(tester)) {}

} // namespace convergence

} // namespace bart
