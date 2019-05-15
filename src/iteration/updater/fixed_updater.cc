#include <data/system.h>
#include "iteration/updater/fixed_updater.h"
#include "formulation/cfem_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template <typename StamperType>
FixedUpdater<StamperType>::FixedUpdater(
    std::unique_ptr<StamperType> stamper_ptr) {
  AssertThrow(stamper_ptr != nullptr,
      dealii::ExcMessage("Error in constructor of FixedUpdater: "
                         "stamper pointer is null."));
  stamper_ptr_ = std::move(stamper_ptr);
}

template <>
void FixedUpdater<formulation::CFEMStamperI>::UpdateFixedTerms(
    data::System& system,
    data::system::GroupNumber group,
    data::system::AngleIndex angle) {

  auto fixed_matrix_ptr_ =
      system.left_hand_side_ptr_->GetFixedTermPtr({group, angle});

  AssertThrow(fixed_matrix_ptr_ != nullptr,
      dealii::ExcMessage("Error in FixedUpdater::UpdateFixedTerms, "
                         "fixed term pointer from lhs is null"));

  *fixed_matrix_ptr_ = 0;

  stamper_ptr_->StampStreamingTerm(*fixed_matrix_ptr_, group);
  stamper_ptr_->StampCollisionTerm(*fixed_matrix_ptr_, group);
  stamper_ptr_->StampBoundaryTerm(*fixed_matrix_ptr_);
}

template class FixedUpdater<formulation::CFEMStamperI>;

} // namespace updater

} // namespace iteration

} // namespace bart