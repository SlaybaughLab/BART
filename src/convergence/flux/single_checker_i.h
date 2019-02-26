#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_

#include "../../data/vector_parameters.h"

namespace bart {

namespace convergence {
 
namespace flux {

/*! \brief Checks for convergence between two provided fluxes. */

class SingleCheckerI {
 public:
  virtual ~SingleCheckerI() = default;
  /* \brief Checks for convergence of two provided fluxes */
  virtual bool CheckIfConverged(data::Flux &, data::Flux &) = 0;
  /* \brief Returns status of convergence (from last call to CheckIfConverged */
  virtual bool is_converged() const = 0;
};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_
