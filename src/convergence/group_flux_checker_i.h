#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_I_H_

#include "../data/vector_parameters.h"

namespace bart {

namespace convergence {

/*! \brief Checks that all group fluxes have converged from the previous
 * iteration. */

class GroupFluxCheckerI {
 public:
  virtual ~GroupFluxCheckerI() = default;
  virtual bool CheckIfConverged(data::ScalarGroupFluxPtrs &current,
                                data::ScalarGroupFluxPtrs &last) = 0;
  virtual bool CheckIfConverged(data::AngularGroupFluxPtrs &current,
                                data::AngularGroupFluxPtrs &last) = 0;
};
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_I_H_
