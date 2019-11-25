#include "formulation/angular/cfem_self_adjoint_angular_flux.h"

namespace bart {

namespace formulation {

namespace angular {

template class CFEMSelfAdjointAngularFlux<1>;
template class CFEMSelfAdjointAngularFlux<2>;
template class CFEMSelfAdjointAngularFlux<3>;

} // namespace angular

} // namespace formulation

} // namespace bart
