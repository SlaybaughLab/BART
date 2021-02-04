#ifndef BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_I_H_
#define BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_I_H_

#include "system/system.hpp"
#include "system/system_types.h"
#include "quadrature/quadrature_types.h"

namespace bart {

namespace formulation {

namespace updater {

class FixedUpdaterI {
 public:
  virtual ~FixedUpdaterI() = default;
  virtual void UpdateFixedTerms(system::System& to_update,
                                system::EnergyGroup,
                                quadrature::QuadraturePointIndex) = 0;
};

} // namespace updater

} // namespace formulation

} // namespace bart

#endif //BART_SRC_FORMULATION_UPDATER_FIXED_UPDATER_I_H_
