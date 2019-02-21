#ifndef BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_

#include "group_fluxes_i.h"

namespace bart {

namespace convergence {

/*! \brief Checks each flux group sequentially for convergence.
 * Will not re-check converged groups until final convergence is believed to
 * have been reached. */

class GroupFluxesSequential : public GroupFluxesI {
 public:
  GroupFluxesSequential() = default;
  ~GroupFluxesSequential() = default;
  bool isConverged(data::GroupFluxes &) {return false; };
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUXES_SEQUENTIAL_H_
