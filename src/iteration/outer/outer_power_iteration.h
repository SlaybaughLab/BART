#ifndef BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_

#include "formulation/updater/fission_source_updater_i.h"
#include "eigenvalue/k_effective/k_effective_updater_i.h"
#include "iteration/outer/outer_iteration.h"

namespace bart {

namespace iteration {

namespace outer {

class OuterPowerIteration : public OuterIteration<double> {
 public:
  using ConvergenceChecker = convergence::FinalI<double>;
  using K_EffectiveUpdater = eigenvalue::k_effective::K_EffectiveUpdaterI;
  using SourceUpdaterType = formulation::updater::FissionSourceUpdaterI;
  using OuterIteration<double>::Reporter;

  OuterPowerIteration(
      std::unique_ptr<GroupIterator> group_iterator_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      std::unique_ptr<K_EffectiveUpdater> k_effective_updater_ptr,
      const std::shared_ptr<SourceUpdaterType> &source_updater_ptr,
      const std::shared_ptr<Reporter> &reporter_ptr = nullptr);
  virtual ~OuterPowerIteration() = default;

  SourceUpdaterType* source_updater_ptr() const {
    return source_updater_ptr_.get();
  };

  K_EffectiveUpdater* k_effective_updater_ptr() const {
    return k_effective_updater_ptr_.get();
  }

 protected:
  convergence::Status CheckConvergence(system::System &system) override;
  void UpdateSystem(system::System &system, const int group, const int angle) override;

  std::shared_ptr<SourceUpdaterType> source_updater_ptr_ = nullptr;
  std::unique_ptr<K_EffectiveUpdater> k_effective_updater_ptr_ = nullptr;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_POWER_ITERATION_H_
