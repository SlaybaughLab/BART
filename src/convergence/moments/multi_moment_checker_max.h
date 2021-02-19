#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_

#include <memory>

#include "convergence/moments/multi_moment_checker.h"

namespace bart {

namespace convergence {

namespace moments {

class MultiMomentCheckerMax : public MultiMomentChecker {
 public:
  using MultiMomentChecker::SingleMomentChecker;
  MultiMomentCheckerMax(std::unique_ptr<SingleMomentChecker> checker) : MultiMomentChecker(std::move(checker)) {}

  bool IsConverged(const system::moments::MomentsMap &current_iteration,
                   const system::moments::MomentsMap &previous_iteration) override;
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_