#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_

#include <memory>

#include "flux_checker_i.h"
#include "group_flux_checker_i.h"


namespace bart {

namespace convergence {

/*! \brief Checks each flux group sequentially for convergence.
 * Will not re-check converged groups until final convergence is believed to
 * have been reached. */

class GroupFluxCheckerSequential : public GroupFluxCheckerI {
 public:
  GroupFluxCheckerSequential(std::unique_ptr<FluxCheckerI> &tester);
  ~GroupFluxCheckerSequential() = default;
  bool isConverged(data::GroupFluxes &current, data::GroupFluxes &last) {
    return false; };

 private:
  /*! Flux convergence tester that will be used to check sequentially */
  std::unique_ptr<FluxCheckerI> tester_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_
