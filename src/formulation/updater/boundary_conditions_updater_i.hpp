#ifndef BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_HPP_
#define BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_HPP_

#include "system/system.hpp"
#include "system/system_types.h"
#include "quadrature/quadrature_types.h"
#include "utility/has_value.hpp"
#include "instrumentation/port.hpp"

namespace bart::formulation::updater {

namespace data_port {
using AggregatedBoundaryConditionValue = instrumentation::Port<double, struct BoundaryConditionParameter>;
} // namespace data_port

class BoundaryConditionsUpdaterI : public utility::HasValue<double>,
                                   public data_port::AggregatedBoundaryConditionValue {
 public:
  virtual ~BoundaryConditionsUpdaterI() = default;
  virtual void UpdateBoundaryConditions(system::System& to_update,
                                        system::EnergyGroup,
                                        quadrature::QuadraturePointIndex) = 0;
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_BOUNDARY_CONDITIONS_UPDATER_I_HPP_
