#ifndef BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_
#define BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_

#include "checker_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks for convergence between two fluxes using some threshold value */
class CheckerThresholdI : public CheckerI {
 public:
  virtual ~CheckerThresholdI() = default;
  /* \brief Set the threshold value for convergence check */
  virtual double SetThreshold() = 0;
  /* \brief Get the threshold value for convergence check */
  virtual double threshold() const = 0;  
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif  // BART_SRC_CONVERGENCE_FLUX_CHECKER_THRESHOLD_I_H_
