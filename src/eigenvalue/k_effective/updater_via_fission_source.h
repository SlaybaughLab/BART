#ifndef BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_H_
#define BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_H_

#include <memory>
#include "calculator/cell/total_aggregated_fission_source_i.h"
#include "eigenvalue/k_effective/updater_via_fission_source_i.h"

namespace bart {

namespace eigenvalue {

namespace k_effective {

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
 * implementation uses a calculator::cell::TotalAggreatedFissionSourceI to
 * accomplish this. With each call of Calculate, the current and previous
 * fission source values are updated. Previous values are not stored in the
 * default implementation.
 *
 */
class UpdaterViaFissionSource : public UpdaterViaFissionSourceI {
 public:
  using FissionSourceCalculator = calculator::cell::TotalAggregatedFissionSourceI;

  UpdaterViaFissionSource(
      std::unique_ptr<FissionSourceCalculator> fission_source_calculator)
      : fission_source_calculator_(std::move(fission_source_calculator)) {}

  /*! \brief Returns the last calculated k_effective */
  double k_effective() const override { return k_effective_; }
  /*! \brief Returns the fission source used in the numerator of the
   * calculation.  */
  double current_fission_source() const { return current_fission_source_; }
  /*! \brief Returns the fission source used in the denominator of the
   * calculation.  */
  double previous_fission_source() const { return previous_fission_source_; }

  /*! \brief Returns raw pointer to dependent fission source calculator */
  FissionSourceCalculator* fission_source_calculator() const {
    return fission_source_calculator_.get();
  }
 private:
  std::unique_ptr<FissionSourceCalculator> fission_source_calculator_;
  double k_effective_ = 1.0;
  double current_fission_source_ = 0;
  double previous_fission_source_ = 0;

};

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_H_