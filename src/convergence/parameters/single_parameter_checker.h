#ifndef BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_
#define BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_

#include <cstdlib>

#include "convergence/single_checker.h"

namespace bart {

namespace convergence {

namespace parameters {

/*! \brief Checks for convergence of any parameter expressed as a double.
 *
 * Convergence of a value \f$x\f$ after the \f$i\f$th iteration is determined by
 * a percentage change from the current iteration value:
 *
 * \f[
 * \Delta_i = \frac{|x_i - x_{i-1}|}{|x_i|}
 * \f]
 *
 * Convergence is achieved if \f$ \Delta_i \leq \Delta_{\text{max}}\f$.
 *
 */

class SingleParameterChecker : public SingleChecker<double> {
 public:
  explicit SingleParameterChecker(double max_delta = 1e-6) {
    max_delta_ = max_delta;
  }

  bool CheckIfConverged(const double &current_iteration,
                        const double &previous_iteration) override {
    double diff = std::abs(current_iteration - previous_iteration);
    return (diff/std::abs(current_iteration) <= max_delta_);
  }

 protected:
  using SingleChecker<double>::max_delta_;
};

} // namespace parameters

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_