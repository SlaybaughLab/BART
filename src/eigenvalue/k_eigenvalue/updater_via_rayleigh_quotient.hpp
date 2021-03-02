#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_

#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

namespace bart::eigenvalue::k_eigenvalue {
/*! \brief Calculates the _k_-eigenvalue using the Rayleigh quotient.
 *
 * The Rayleigh quotient provides an updated value of k:
 * \f[
 * k_{i +1} = k_{i}\frac{\phi_{i}^T\phi_{i + 1}}{\phi_{i}^T\phi_{i}}
 * \f]
 *
 */
class UpdaterViaRayleighQuotient : public K_EigenvalueCalculatorI {
 public:
  UpdaterViaRayleighQuotient() { this->set_description("k-effective updater via Rayleigh quotient. "); };
  [[nodiscard]] auto CalculateK_Eigenvalue(system::System &system) -> double override;
  [[nodiscard]] auto k_eigenvalue() const -> std::optional<double> override { return last_calculated_k_effective_; };
 private:
  std::optional<double> last_calculated_k_effective_{ std::nullopt };
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif //BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_
