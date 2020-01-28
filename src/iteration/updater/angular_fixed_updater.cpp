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

template class AngularFixedUpdater<formulation::AngularStamperI<1>>;
template class AngularFixedUpdater<formulation::AngularStamperI<2>>;
template class AngularFixedUpdater<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
