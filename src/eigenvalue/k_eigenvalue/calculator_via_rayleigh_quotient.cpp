#include "eigenvalue/k_eigenvalue/calculator_via_rayleigh_quotient.hpp"

namespace bart::eigenvalue::k_eigenvalue {

double CalculatorViaRayleighQuotient::CalculateK_Eigenvalue(system::System& system) {
  const double current_k_effective{ system.k_effective.value_or(1.0) };
  double calculated_k_effective{ current_k_effective };

  dealii::Vector<double> current_flux = system.current_moments->GetMoment(std::array{0, 0, 0});
  dealii::Vector<double> previous_flux = system.previous_moments->GetMoment(std::array{0, 0, 0});

  for (int group = 1; group < system.total_groups; ++group) {
    std::array<int, 3> index{ group, 0, 0 };
    current_flux.add(1, system.current_moments->GetMoment(index));
    previous_flux.add(1, system.previous_moments->GetMoment(index));
  }

  if (!previous_flux.all_zero()) {
    calculated_k_effective = current_k_effective * (current_flux * previous_flux) / (previous_flux * previous_flux);
  }

  last_calculated_k_effective_ = calculated_k_effective;
  return calculated_k_effective;
}

} // namespace bart::eigenvalue::k_eigenvalue
