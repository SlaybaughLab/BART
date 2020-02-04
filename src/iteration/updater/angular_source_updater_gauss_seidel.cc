#include "angular_source_updater_gauss_seidel.h"

#include "formulation/angular_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template<typename StamperType>
AngularSourceUpdaterGaussSeidel<StamperType>::AngularSourceUpdaterGaussSeidel(
    std::shared_ptr<StamperType> stamper_ptr,
    std::shared_ptr<QuadratureSetType> quadrature_set_ptr)
    : SourceUpdater<StamperType>(stamper_ptr) {
  AssertThrow(stamper_ptr != nullptr,
      dealii::ExcMessage("Error in constructor of AngularSourceUpdaterGaussSeidel,"
                         " stamper pointer is null"))
  AssertThrow(quadrature_set_ptr != nullptr,
              dealii::ExcMessage("Error in constructor of AngularSourceUpdaterGaussSeidel,"
                                 " quadrature set pointer is null"))
}

template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<1>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<2>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
