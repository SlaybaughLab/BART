#ifndef BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_

#include <memory>

#include "flux_i.h"
#include "group_fluxes_i.h"


namespace bart {

namespace convergence {

/*! \brief Checks each flux group sequentially for convergence.
 * Will not re-check converged groups until final convergence is believed to
 * have been reached. */

class GroupFluxesSequential : public GroupFluxesI {
 public:
  GroupFluxesSequential(std::unique_ptr<FluxI> &tester);
  ~GroupFluxesSequential() = default;
  bool isConverged(data::GroupFluxes &) {return false; };

 private:
  /*! Flux convergence tester that will be used to check sequentially */
  std::unique_ptr<FluxI> tester_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_
