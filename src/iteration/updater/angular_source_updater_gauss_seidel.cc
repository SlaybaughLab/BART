#include "angular_source_updater_gauss_seidel.h"

#include "formulation/angular_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<1>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<2>>;
template class AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
