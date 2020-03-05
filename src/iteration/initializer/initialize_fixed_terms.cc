#include "iteration/initializer/initialize_fixed_terms.h"

namespace bart {

namespace iteration {

namespace initializer {

InitializeFixedTerms::InitializeFixedTerms(
    std::unique_ptr<FixedUpdaterType> fixed_updater_ptr,
    int total_groups,
    int total_angles)
    : fixed_updater_ptr_(std::move(fixed_updater_ptr)),
      total_groups_(total_groups),
      total_angles_(total_angles) {
  AssertThrow(total_groups_ > 0,
              dealii::ExcMessage("Error in iteration::initializer::InitializeFixedTerms "
                                 "constructor: total_groups !> 0"));
  AssertThrow(total_angles_ > 0,
              dealii::ExcMessage("Error in iteration::initializer::InitializeFixedTerms "
                                 "constructor: total_angles !> 0"));
  AssertThrow(fixed_updater_ptr_ != nullptr,
              dealii::ExcMessage("Error in iteration::initializer::InitializeFixedTerms "
                                 "constructor: fixed updater ptr is nullptr"))
}

void InitializeFixedTerms::Initialize(system::System &sys) {
  for (int group = 0; group < total_groups_; ++group) {
    for (int angle = 0; angle < total_angles_; ++angle) {
      system::EnergyGroup energy_group(group);
      quadrature::QuadraturePointIndex angle_index(angle);
      fixed_updater_ptr_->UpdateFixedTerms(sys, energy_group, angle_index);
    }
  }
}

} // namespace initializer

} // namespace iteration

} // namespace bart
