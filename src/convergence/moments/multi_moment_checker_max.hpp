#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_

#include <memory>

#include "convergence/moments/multi_moment_checker.hpp"

namespace bart::convergence::moments {

class MultiMomentCheckerMax : public MultiMomentChecker {
 public:
  using MultiMomentChecker::SingleMomentChecker;
  using MomentsMap = system::moments::MomentsMap;
  MultiMomentCheckerMax(std::unique_ptr<SingleMomentChecker> checker) : MultiMomentChecker(std::move(checker)) {}

  bool IsConverged(const MomentsMap &current_iteration, const MomentsMap &previous_iteration) override;
};

} // namespace bart::convergence::moments

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_