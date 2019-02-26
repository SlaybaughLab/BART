#ifndef BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_
#define BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_

#include "single_checker_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks for convergence between two fluxes using the value of their
 * difference */
class SingleCheckerDeltaI : public SingleCheckerI {
 public:
  virtual ~SingleCheckerDeltaI() = default;
  /* \brief Set the threshold value for convergence check */
  virtual double SetMaxDelta() = 0;
  /* \brief Get the threshold value for convergence check */
  virtual double max_delta() const = 0;
  /* \brief Get the delta value from the last convergence check */
  virtual double delta() const = 0;
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif  // BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_
