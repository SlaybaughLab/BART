#ifndef BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_HPP_

#include <cstdlib>

#include "convergence/convergence_checker.hpp"

namespace bart::convergence::parameters {

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

class SingleParameterChecker : public ConvergenceChecker<double, double> {
 public:
  explicit SingleParameterChecker(double max_delta = 1e-6) { SetMaxDelta(max_delta); }

  auto SetMaxDelta(const double& to_set) -> void override {
    AssertThrow(to_set > 0, dealii::ExcMessage("Error in SingleParameterChecker::SetMaxDelta, value to set must be "
                                               "greater than 0"))
    ConvergenceChecker<double, double>::SetMaxDelta(to_set);
  }

  [[nodiscard]] auto CheckIfConverged(const double &current_value, const double &previous_value) -> bool override {
    double diff = std::abs(current_value - previous_value);
    this->delta_ = diff;
    return (diff/std::abs(current_value) <= max_delta_);
  }

 protected:
  using ConvergenceChecker<double>::max_delta_;
};

} // namespace bart::convergence::parameters

#endif // BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_HPP_