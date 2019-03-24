#ifndef BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_
#define BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_

#include "convergence/single_checker.h"

namespace bart {

namespace convergence {

namespace parameters {

class SingleParameterChecker : public SingleChecker<double> {
 public:
  explicit SingleParameterChecker(double max_delta = 1e-6) {
    max_delta_ = max_delta;
  }

  bool CheckIfConverged(const double &current_iteration,
                        const double &previous_iteration) override {
    return false;
  }

 protected:
  using SingleChecker<double>::max_delta_;
};

} // namespace parameters

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_PARAMETERS_SINGLE_PARAMETER_CHECKER_H_