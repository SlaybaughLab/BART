#include "iteration/outer/outer_power_iteration.hpp"

namespace bart {

namespace iteration {

namespace outer {

OuterPowerIteration::OuterPowerIteration(
    std::unique_ptr<GroupIterator> group_iterator_ptr,
    std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
    std::unique_ptr<K_EffectiveUpdater> k_effective_updater_ptr,
    const std::shared_ptr<SourceUpdaterType> &source_updater_ptr)
    : OuterIteration(
        std::move(group_iterator_ptr),
        std::move(convergence_checker_ptr)),
      source_updater_ptr_(source_updater_ptr),
      k_effective_updater_ptr_(std::move(k_effective_updater_ptr)) {

  AssertThrow(k_effective_updater_ptr_ != nullptr,
              dealii::ExcMessage("KEffective updater pointer passed to "
                                 "OuterPowerIteration constructor is null"));
  AssertThrow(source_updater_ptr_ != nullptr,
              dealii::ExcMessage("Source updater pointer passed to OuterIteration "
                                 "constructor is null"));
  this->set_description("outer power iteration",
                        utility::DefaultImplementation(true));
}

convergence::Status OuterPowerIteration::CheckConvergence(system::System &system) {

  double k_effective_last = system.k_effective.value_or(0.0);
  system.k_effective = k_effective_updater_ptr_->CalculateK_Eigenvalue(system);

  return convergence_checker_ptr_->ConvergenceStatus(
      system.k_effective.value(), k_effective_last);
}
void OuterPowerIteration::UpdateSystem(system::System &system, const int group,
    const int angle) {
  source_updater_ptr_->UpdateFissionSource(system,
      system::EnergyGroup(group),
      quadrature::QuadraturePointIndex(angle));
}

} // namespace outer

} // namespace iteration

} // namespace bart