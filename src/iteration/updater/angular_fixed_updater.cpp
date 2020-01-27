#include "iteration/updater/angular_fixed_updater.h"
#include "formulation/cfem_saaf_stamper.h"

namespace bart {

namespace iteration {

namespace updater {

template class AngularFixedUpdater<formulation::CFEM_SAAF_Stamper<1>>;
template class AngularFixedUpdater<formulation::CFEM_SAAF_Stamper<2>>;
template class AngularFixedUpdater<formulation::CFEM_SAAF_Stamper<3>>;

} // namespace updater

} // namespace iteration

} // namespace bart
