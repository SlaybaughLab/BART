#include "iteration/updater/angular_fixed_updater.h"
#include "formulation/angular_stamper_i.h"

namespace bart {

namespace iteration {

namespace updater {

template class AngularFixedUpdater<formulation::AngularStamperI<1>>;
template class AngularFixedUpdater<formulation::AngularStamperI<2>>;
template class AngularFixedUpdater<formulation::AngularStamperI<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
