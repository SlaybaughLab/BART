#include "iteration/initializer/initialize_fixed_terms.h"

namespace bart {

namespace iteration {

namespace initializer {

InitializeFixedTerms::InitializeFixedTerms(
    std::unique_ptr<FixedUpdaterType> fixed_updater_ptr)
    : fixed_updater_ptr_(std::move(fixed_updater_ptr)) {
  AssertThrow(fixed_updater_ptr_ != nullptr,
              dealii::ExcMessage("Error in iteration::initializer::InitializeFixedTerms "
                                 "constructor: fixed updater ptr is nullptr"))
}

void InitializeFixedTerms::Initialize(system::System &sys) {

}

} // namespace initializer

} // namespace iteration

} // namespace bart
