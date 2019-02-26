#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_H_

#include "single_checker_delta_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks for convergence between two fluxes using the value of their
 * difference */
class SingleCheckerDelta : public SingleCheckerDeltaI {
 public:
  virtual ~SingleCheckerDelta() = default;
  void SetMaxDelta(double to_set) override { max_delta_ = to_set; }
  double max_delta() const override { return max_delta_; }
  double delta() const override { return delta_; }
 protected:
  double max_delta_ = 1e-6;
  double delta_ = 0;
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif  // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_H_
