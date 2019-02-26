#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_

#include "single_checker_i.h"
#include "../../data/vector_parameters.h"

namespace bart {

namespace convergence {
 
namespace flux {

/*! \brief Checks for convergence between two provided fluxes. */

class SingleChecker : public SingleCheckerI {
 public:
  virtual ~SingleChecker() = default;
  bool is_converged() const override { return is_converged_; };
 protected:
  bool is_converged_ = false;
};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_
