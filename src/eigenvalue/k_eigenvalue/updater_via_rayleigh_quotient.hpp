#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_

#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

namespace bart::eigenvalue::k_eigenvalue {

class UpdaterViaRayleighQuotient : public K_EigenvalueCalculatorI {
 public:
  UpdaterViaRayleighQuotient() { this->set_description("k-effective updater via Rayleigh quotient. "); };
  double CalculateK_Eigenvalue(system::System &system) override;
  std::optional<double> k_eigenvalue() const override {
    return last_calculated_k_effective_; };
 private:
  std::optional<double> last_calculated_k_effective_{ std::nullopt };
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif //BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_RAYLEIGH_QUOTIENT_HPP_
