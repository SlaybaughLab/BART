#ifndef BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_H_
#define BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_H_

#include "system/system.h"
#include "system/system_types.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace formulation {

namespace updater {

class BoundaryConditionsUpdaterI {
 public:
  virtual ~BoundaryConditionsUpdaterI() = default;
  virtual void UpdateBoundaryConditions(system::System& to_update,
                                        system::EnergyGroup,
                                        quadrature::QuadraturePointIndex) = 0;
};


} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_H_
