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

}

} // namespace initializer

} // namespace iteration

} // namespace bart
