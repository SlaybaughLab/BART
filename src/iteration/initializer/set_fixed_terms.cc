#include "iteration/initializer/set_fixed_terms.h"

#include "system/system.h"
#include "iteration/updater/fixed_updater_i.h"

namespace bart {

namespace iteration {

namespace initializer {

SetFixedTerms::SetFixedTerms(
    std::unique_ptr<updater::FixedUpdaterI> fixed_updater_ptr,
    const int total_groups,
    const int total_angles)
    : fixed_updater_ptr_(std::move(fixed_updater_ptr)),
      total_groups_(total_groups),
      total_angles_(total_angles) {

  AssertThrow(total_groups_ > 0,
              dealii::ExcMessage("Error in iteration::initializer::SetFixedTerms "
                                 "constructor: total_groups !> 0"));
  AssertThrow(total_angles_ > 0,
              dealii::ExcMessage("Error in iteration::initializer::SetFixedTerms "
                                 "constructor: total_angles !> 0"));
  AssertThrow(fixed_updater_ptr_ != nullptr,
              dealii::ExcMessage("Error in iteration::initializer::SetFixedTerms "
                                 "constructor: fixed updater ptr is nullptr"));
}

void SetFixedTerms::Initialize(system::System &sys) {
  for (int group = 0; group < total_groups_; ++group) {
    for (int angle = 0; angle < total_angles_; ++angle) {
      fixed_updater_ptr_->UpdateFixedTerms(sys, group, angle);
    }
  }
}

} // namespace initializer

} // namespace iteration

} // namespace bart