#ifndef BART_SRC_FORMULATION_UPDATER_SCATTERING_SOURCE_UPDATER_I_H_
#define BART_SRC_FORMULATION_UPDATER_SCATTERING_SOURCE_UPDATER_I_H_

#include "system/system.hpp"
#include "system/system_types.h"
#include "quadrature/quadrature_types.h"
#include "utility/has_value.hpp"
#include "instrumentation/port.hpp"

namespace bart::formulation::updater {

namespace data_port {
using AggregatedScatteringSourceValue = instrumentation::Port<double, struct ScatteringSourceValueParameter>;
} // namespace data_port

 class ScatteringSourceUpdaterI : public utility::HasValue<double>, public data_port::AggregatedScatteringSourceValue {
 public:
  virtual ~ScatteringSourceUpdaterI() = default;
  virtual void UpdateScatteringSource(system::System& to_update,
                                      system::EnergyGroup,
                                      quadrature::QuadraturePointIndex) = 0;
};

} // namespace bart::formulation::updater

#endif //BART_SRC_FORMULATION_UPDATER_SCATTERING_SOURCE_UPDATER_I_H_
