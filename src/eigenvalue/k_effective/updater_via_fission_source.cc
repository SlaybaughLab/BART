#include "eigenvalue/k_effective/updater_via_fission_source.h"

namespace bart {

namespace eigenvalue {

namespace k_effective {

UpdaterViaFissionSource::UpdaterViaFissionSource(
    std::unique_ptr<FissionSourceCalculator> fission_source_calculator,
    double initial_k_effective,
    double initial_fission_source)
    : fission_source_calculator_(std::move(fission_source_calculator)),
      initial_k_effective_(initial_k_effective),
      initial_fission_source_(initial_fission_source) {}

double UpdaterViaFissionSource::CalculateK_Effective(system::System &system) {

  current_fission_source_ =
      fission_source_calculator_->AggregatedFissionSource(
          system.current_moments.get());

  k_effective_ = initial_k_effective_ *
      current_fission_source_.value() / initial_fission_source_;

  return k_effective_.value();
}

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart