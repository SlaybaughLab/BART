#ifndef BART_SRC_CONVERGENCE_MOMENTS_CONVERGENCE_CHECKER_L_INFINITY_NORM_HPP_
#define BART_SRC_CONVERGENCE_MOMENTS_CONVERGENCE_CHECKER_L_INFINITY_NORM_HPP_

#include "convergence/convergence_checker.hpp"
#include <deal.II/lac/vector.h>

namespace bart::convergence::moments {

class ConvergenceCheckerLInfinityNorm : public ConvergenceChecker<dealii::Vector<double>> {
 public:
  using Vector = dealii::Vector<double>;
  explicit ConvergenceCheckerLInfinityNorm(const double max_delta = 1e-6) { max_delta_ = CheckNonNegative(max_delta); };

  auto SetMaxDelta(const double& to_set) -> void override { max_delta_ = CheckNonNegative(to_set); };
  auto IsConverged(const Vector &current_iteration, const Vector &previous_iteration) -> bool override {
    Vector difference(current_iteration);
    difference.add(-1, previous_iteration);
    delta_ = difference.linfty_norm()/current_iteration.linfty_norm();
    is_converged_ = delta_ <= max_delta_;
    return is_converged_;
  }

 private:
  auto CheckNonNegative(const double to_check) -> double {
    AssertThrow(to_check > 0, dealii::ExcMessage("Error in ConvergenceCheckerLInfinityNorm, max delta value to set must "
                                                 "be greater than 0"))
    return to_check;
  }
};

} // namespace bart::convergence::moments

#endif //BART_SRC_CONVERGENCE_MOMENTS_CONVERGENCE_CHECKER_L_INFINITY_NORM_HPP_
