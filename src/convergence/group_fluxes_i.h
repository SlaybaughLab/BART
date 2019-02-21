#ifndef BART_SRC_CONVERGENCE_GROUP_FLUXES_I_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUXES_I_H_

#include "../data/vector_parameters.h"

namespace bart {

namespace convergence {

/*! \brief Checks that all fluxes have converged across groups */

class GroupFluxesI {
 public:
  virtual ~GroupFluxesI() = default;
  virtual bool isConverged(data::GroupFluxes &) = 0;
};
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUXES_I_H_
