#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_

#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

namespace bart {

namespace eigenvalue {

namespace k_eigenvalue {

class UpdaterViaFissionSourceI : public K_EigenvalueCalculatorI {
 public:
  virtual ~UpdaterViaFissionSourceI() = default;

  /*! \brief Return the current problem fission source. */
  virtual std::optional<double> current_fission_source() const = 0;
};

} // namespace k_eigenvalue

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_