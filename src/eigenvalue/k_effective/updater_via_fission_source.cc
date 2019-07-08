#include "eigenvalue/k_effective/updater_via_fission_source.h"

namespace bart {

namespace eigenvalue {

namespace k_effective {

double UpdaterViaFissionSource::CalculateK_Effective(system::System &system) {

  previous_fission_source_ = current_fission_source_;

  current_fission_source_ =
      fission_source_calculator_->AggregatedFissionSource(
          system.current_moments.get());

  k_effective_ = system.k_effective.value_or(1.0) *
      current_fission_source_.value() /
      previous_fission_source_.value_or(1.0);

  return k_effective_.value();
}

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart