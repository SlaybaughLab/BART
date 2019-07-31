#ifndef BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_

#include "eigenvalue/k_effective/k_effective_updater_i.h"
#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterPowerIteration : public OuterIteration<double> {
 public:
  using ConvergenceChecker = convergence::FinalI<double>;
  using K_EffectiveUpdater = eigenvalue::k_effective::K_EffectiveUpdaterI;
  using OuterIteration<double>::SourceUpdater;

  OuterPowerIteration(
      std::unique_ptr<GroupIterator> group_iterator_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      std::unique_ptr<K_EffectiveUpdater> k_effective_updater_ptr,
      const std::shared_ptr<SourceUpdater> &source_updater_ptr);
  virtual ~OuterPowerIteration() = default;

  K_EffectiveUpdater* k_effective_updater_ptr() const {
    return k_effective_updater_ptr_.get();
  }

 protected:
  convergence::Status CheckConvergence(system::System &system) override;
  void UpdateSystem(system::System &system, const int group, const int angle) override;

  std::unique_ptr<K_EffectiveUpdater> k_effective_updater_ptr_ = nullptr;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
