#include "iteration/updater/angular_fixed_updater.h"
#include "formulation/angular_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template<typename StamperType>
AngularFixedUpdater<StamperType>::AngularFixedUpdater(
    std::shared_ptr<StamperType> stamper_ptr,
    std::shared_ptr<QuadratureSetType> quadrature_set_ptr)
    : stamper_ptr_(stamper_ptr),
      quadrature_set_ptr_(quadrature_set_ptr) {
  AssertThrow(stamper_ptr != nullptr,
      dealii::ExcMessage("Error in constructor of AngularFixedUpdater, "
                         "stamper pointer is null"));
  AssertThrow(quadrature_set_ptr != nullptr,
              dealii::ExcMessage("Error in constructor of AngularFixedUpdater, "
                                 "quadrature set pointer is null"));
}
template<typename StamperType>
void AngularFixedUpdater<StamperType>::UpdateFixedTerms(
    system::System &system,
    system::GroupNumber group,
    system::AngleIndex angle) {

  auto fixed_matrix_ptr_ =
      system.left_hand_side_ptr_->GetFixedTermPtr({group, angle});

  AssertThrow(fixed_matrix_ptr_ != nullptr,
              dealii::ExcMessage("Error in AngularFixedUpdater::UpdateFixedTerms"
                                 ", fixed term pointer from lhs is null"));

  *fixed_matrix_ptr_ = 0;

  auto quadrature_point_ptr = quadrature_set_ptr_->GetQuadraturePoint(
      quadrature::QuadraturePointIndex(angle));

  stamper_ptr_->StampStreamingTerm(*fixed_matrix_ptr_,
                                   quadrature_point_ptr,
                                   system::EnergyGroup(group));
  stamper_ptr_->StampCollisionTerm(*fixed_matrix_ptr_,
                                   system::EnergyGroup(group));
  stamper_ptr_->StampBoundaryBilinearTerm(*fixed_matrix_ptr_,
                                          quadrature_point_ptr,
                                          system::EnergyGroup(group));
}

template class AngularFixedUpdater<formulation::AngularStamperI<1>>;
template class AngularFixedUpdater<formulation::AngularStamperI<2>>;
template class AngularFixedUpdater<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
