#ifndef BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_
#define BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_

#include "convergence/moments/single_moment_checker_i.h"

namespace bart {

namespace convergence {

namespace moments {

/*! \brief Checks for convergence between flux moments using the percentage
 * change in the L1 norm, compared to the current iteration (\f$i\f$):
 *
 * \f[
 *
 * \Delta_i = \frac{|\phi_i - \phi_{i-1}|_{1}}{|\phi_{i}|_{1}}
 *
 * \f]
 *
 * Convergence is achieved if \f$\Delta_i \leq \Delta_{\text{max}}\f$.
 * */

class SingleMomentCheckerL1Norm : public SingleMomentCheckerI {
 public:
  /*! \brief Default constructor, setting max delta to \f$10^{-6}\f$. */

  explicit SingleMomentCheckerL1Norm(const double max_delta = 1e-6) {
    max_delta_ = max_delta;
  };

  auto SetMaxDelta(const double& to_set) -> void override {
    AssertThrow(to_set > 0, dealii::ExcMessage("Error in SingleParameterChecker::SetMaxDelta, value to set must be "
                                               "greater than 0"))
    SingleMomentCheckerI::SetMaxDelta(to_set);
  }

  ~SingleMomentCheckerL1Norm() = default;

  bool CheckIfConverged(
      const system::moments::MomentVector &current_iteration,
      const system::moments::MomentVector &previous_iteration) override;
};

} // namespace moments
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_
