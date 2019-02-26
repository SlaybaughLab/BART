#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_

#include "single_checker_delta_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks for convergence between two fluxes using the percentage
 * change in the L1 norms */

class SingleCheckerDeltaL1Norm : public SingleCheckerDeltaI {
 public:
  SingleCheckerDeltaL1Norm() = default;
  ~SingleCheckerDeltaL1Norm() = default;
  bool CheckIfConverged(data::Flux &current_iteration,
                        data::Flux &previous_iteration) override;
  bool is_converged() const { return is_converged_; } override;
  void SetMaxDelta(double to_set) { max_delta_ = to_set; } override;
  double max_delta() const { return max_delta_; };
  double delta() const { return delta_; };  

 private:
  bool is_converged_ = false;;
  double max_delta_ = 1e-6;
  double delta_ = 1e-3;
};
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_DELTA_L1_NORM_H_
