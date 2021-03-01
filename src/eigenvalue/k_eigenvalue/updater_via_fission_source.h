#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_H_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_H_

#include <memory>
#include <optional>
#include "calculator/cell/total_aggregated_fission_source_i.hpp"
#include "eigenvalue/k_eigenvalue/updater_via_fission_source_i.h"

namespace bart {

namespace eigenvalue {

namespace k_eigenvalue {

/*! \brief Calculates an updated k eigenvalue using fission source.
 *
 * This class calculates an updated k eigenvalue using the following equation,
 * \f[
 *
 * k_{\text{eff}}^{k+1} = k_{\text{eff}}^k\frac{\int\nu\Sigma_f\phi^{k+1}}{\int\nu\Sigma_f\phi^{k}}\;,
 *
 * \f]
 *
 * where the integration is performed over the entire neutron phase space. This
 * implementation uses a calculator::cell::TotalAggregatedFissionSourceI to
 * accomplish this. With each call of CalculateKEff, the current and previous
 * fission source values are updated. Previous values are not stored in the
 * default implementation.
 *
 */
class UpdaterViaFissionSource : public UpdaterViaFissionSourceI {
 public:
  using FissionSourceCalculator = calculator::cell::TotalAggregatedFissionSourceI;

  UpdaterViaFissionSource(
      std::unique_ptr<FissionSourceCalculator> fission_source_calculator,
      double initial_k_effective,
      double initial_fission_source);

  double CalculateK_Eigenvalue(system::System& system) override;
  /*! \brief Returns the last calculated k_effective */
  std::optional<double> k_eigenvalue() const override {
    return k_effective_; }

  /*! \brief Returns initial k_effective guess */
  double initial_k_eigenvalue() const {
    return initial_k_effective_;
  }

  /*! \brief Returns the fission source used in the numerator of the
   * calculation.  */
  std::optional<double> current_fission_source() const override {
    return current_fission_source_; }
  /*! \brief Returns the fission source used in the denominator of the
   * calculation.  */
  std::optional<double> initial_fission_source() const {
    return initial_fission_source_; }

  /*! \brief Returns raw pointer to dependent fission source calculator */
  FissionSourceCalculator* fission_source_calculator() const {
    return fission_source_calculator_.get();
  }
 private:
  std::unique_ptr<FissionSourceCalculator> fission_source_calculator_;
  std::optional<double> k_effective_ = std::nullopt;
  std::optional<double> current_fission_source_ = std::nullopt;
  const double initial_k_effective_;
  const double initial_fission_source_;

};

} // namespace k_eigenvalue

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_H_