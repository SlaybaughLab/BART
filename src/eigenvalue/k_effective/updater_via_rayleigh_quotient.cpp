#include "eigenvalue/k_effective/updater_via_rayleigh_quotient.hpp"

namespace bart::eigenvalue::k_effective {

double UpdaterViaRayleighQuotient::CalculateK_Effective(system::System& system) {
  double calculated_k_effective{ 0 };
  bool any_previous_flux_is_nonzero{ false };
  for (int group = 0; group < system.total_groups; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    auto current_flux = system.current_moments->GetMoment(index);
    auto previous_flux = system.previous_moments->GetMoment(index);
    if (!previous_flux.all_zero()) {
      any_previous_flux_is_nonzero = true;
      calculated_k_effective += system.k_effective.value_or(1) * (current_flux * previous_flux) / (previous_flux * previous_flux);
    }
  }

  if (!any_previous_flux_is_nonzero) {
    calculated_k_effective = system.k_effective.value_or(1.0);
  }

  last_calculated_k_effective_ = calculated_k_effective;
  return calculated_k_effective;
}

} // namespace bart::eigenvalue::k_effective
