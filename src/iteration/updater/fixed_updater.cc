#include <system/system.h>
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

/*! \brief Updates the fixed source term for a CFEM formulation.
 *
 * This implementation for the formulation::CFEMStamperI class first retrieves the
 * appropriate matrix for the provided group and angle from left hand side. Then,
 * the matrix is set to zero, and the Streaming, Collision, and Boundary terms
 * are stamped.
 *
 * @param system system to retrieve the left hand side
 * @param group group for stamping streaming and collision terms and
 *              retrieving the correct left hand side
 * @param angle used to retrieve the correct left hand side, not needed for
 *              stamping.
 */
template <>
void FixedUpdater<formulation::CFEMStamperI>::UpdateFixedTerms(
    system::System& system,
    const data::system::GroupNumber group,
    const data::system::AngleIndex angle) {

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