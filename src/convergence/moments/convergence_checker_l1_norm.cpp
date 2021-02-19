#include "convergence/moments/convergence_checker_l1_norm.hpp"

namespace bart::convergence::moments {

ConvergenceCheckerL1Norm::ConvergenceCheckerL1Norm(const double max_delta) {
  max_delta_ = max_delta;
}

bool ConvergenceCheckerL1Norm::IsConverged(const Vector& current_iteration, const Vector& previous_iteration) {
  Vector difference(current_iteration);
  difference -= previous_iteration;
  delta_ = difference.l1_norm()/current_iteration.l1_norm();
  is_converged_ = delta_ <= max_delta_;
  return is_converged_;
}

auto ConvergenceCheckerL1Norm::SetMaxDelta(const double &to_set) -> void {
  AssertThrow(to_set > 0, dealii::ExcMessage("Error in ConvergenceCheckerL1Norm::SetMaxDelta, value to set must be "
                                             "greater than 0"))
  this->SetMaxDelta(to_set);
}

} // namespace bart::convergence::moments