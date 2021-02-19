#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_

#include <memory>

#include "convergence/moments/multi_moment_checker.hpp"

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

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_HPP_