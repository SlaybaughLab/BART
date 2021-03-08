#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_HPP_

#include <memory>
#include <optional>

#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

namespace bart::eigenvalue::k_eigenvalue {

/*! \brief Calculates an updated _k_-eigenvalue using fission source.
 *
 * This class calculates an updated _k_-eigenvalue using the following equation,
 * \f[
 *
 * k_{\text{eff}}^{k+1} = k_{\text{eff}}^k\frac{\int\nu\Sigma_f\phi^{k+1}}{\int\nu\Sigma_f\phi^{k}}\;,
 *
 * \f]
 *
 * where the integration is performed over the entire neutron phase space. This implementation uses a
 * calculator::cell::TotalAggregatedFissionSourceI to accomplish this. With each call of CalculateKEff, the current and
 * previous fission source values are updated. Previous values are not stored in the default implementation.
 */
class CalculatorViaFissionSource : public K_EigenvalueCalculatorI {
 public:
  using FissionSourceCalculator = calculator::cell::TotalAggregatedFissionSourceI;

  /*! \brief Constructor.
   *
   * @param fission_source_calculator calculator for domain fission source
   * @param initial_k_eigenvalue initial guess for _k_-eigenvalue
   * @param initial_fission_source initial guess for fission source
   */
  CalculatorViaFissionSource(std::unique_ptr<FissionSourceCalculator>, double initial_k_eigenvalue,
                             double initial_fission_source);

  [[nodiscard]] auto CalculateK_Eigenvalue(system::System& system) -> double override;
  /*! \brief Returns the last calculated k_eigenvalue */
  [[nodiscard]] auto k_eigenvalue() const -> std::optional<double> override { return k_eigenvalue_; }
  /*! \brief Returns initial k_eigenvalue guess */
  [[nodiscard]] auto initial_k_eigenvalue() const -> double { return initial_k_eigenvalue_; }
  /*! \brief Returns the fission source used in the numerator of the calculation.  */
  [[nodiscard]] auto current_fission_source() const -> std::optional<double> { return current_fission_source_; }
  /*! \brief Returns the fission source used in the denominator of the calculation.  */
  [[nodiscard]] auto initial_fission_source() const -> std::optional<double> { return initial_fission_source_; }

  /*! \brief Returns raw pointer to dependent fission source calculator */
  [[nodiscard]] auto fission_source_calculator() const { return fission_source_calculator_.get(); }
 protected:
  std::unique_ptr<FissionSourceCalculator> fission_source_calculator_;
  std::optional<double> k_eigenvalue_{ std::nullopt };
  std::optional<double> current_fission_source_{ std::nullopt };
  const double initial_k_eigenvalue_;
  const double initial_fission_source_;
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_HPP_