#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_

#include "eigenvalue/k_eigenvalue/k_effective_updater_i.hpp"

namespace bart {

namespace eigenvalue {

namespace k_effective {

class UpdaterViaFissionSourceI : public K_EffectiveUpdaterI {
 public:
  virtual ~UpdaterViaFissionSourceI() = default;

  /*! \brief Return the current problem fission source. */
  virtual std::optional<double> current_fission_source() const = 0;
};

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_UPDATER_VIA_FISSION_SOURCE_I_H_