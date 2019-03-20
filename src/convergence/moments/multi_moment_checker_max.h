#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_

#include "convergence/moments/multi_moment_checker.h"

namespace bart {

namespace convergence {

namespace moments {

class MultiMomentCheckerMax : public MultiMomentChecker {
 public:
  MultiMomentCheckerMax() = default;
  ~MultiMomentCheckerMax() = default;

  bool CheckIfConverged(const data::MomentsMap &current_iteration,
                        const data::MomentsMap &previous_iteration) override {
    return false;
  }
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_