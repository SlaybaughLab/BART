#ifndef BART_SRC_CONVERGENCE_FLUX_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_FLUX_CHECKER_I_H_

#include "../../data/vector_parameters.h"

namespace bart {

namespace convergence {
 
namespace flux {

/*! \brief Checks for convergence between two provided fluxes. */

class CheckerI {
 public:
  virtual ~CheckerI() = default;
  virtual bool CheckIfConverged(data::Flux &, data::Flux &) = 0;
  virtual bool is_converged() const = 0;
};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_CHECKER_I_H_
