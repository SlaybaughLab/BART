#ifndef BART_SRC_FORMULATION_UPDATER_FISSION_SOURCE_UPDATER_I_HPP_
#define BART_SRC_FORMULATION_UPDATER_FISSION_SOURCE_UPDATER_I_HPP_

#include "system/system.h"
#include "system/system_types.h"
#include "quadrature/quadrature_types.h"
#include "utility/has_value.hpp"
#include "instrumentation/port.hpp"

namespace bart::formulation::updater {

namespace data_port {
using AggregatedFissionSourceValue = instrumentation::Port<double, struct FissionSourceValueParameter>;
} // namespace data_port

class FissionSourceUpdaterI : public utility::HasValue<double>, public data_port::AggregatedFissionSourceValue {
 public:
  virtual ~FissionSourceUpdaterI() = default;
  virtual void UpdateFissionSource(system::System& to_update,
                                   system::EnergyGroup,
                                   quadrature::QuadraturePointIndex) = 0;
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_FISSION_SOURCE_UPDATER_I_HPP_
