#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_

#include <memory>

#include "convergence/moments/single_moment_checker_i.h"
#include "convergence/moments/multi_moment_checker.h"

namespace bart {

namespace convergence {

namespace moments {

class MultiMomentCheckerMax : public MultiMomentChecker {
 public:
  MultiMomentCheckerMax(std::unique_ptr<SingleMomentCheckerI> checker)
      : MultiMomentChecker(std::move(checker)) {}
  ~MultiMomentCheckerMax() = default;

  bool CheckIfConverged(const data::MomentsMap &current_iteration,
                        const data::MomentsMap &previous_iteration) override;
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_MAX_H_