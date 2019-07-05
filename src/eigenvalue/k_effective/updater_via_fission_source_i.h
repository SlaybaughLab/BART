#ifndef BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_I_H_
#define BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_I_H_

#include "eigenvalue/k_effective/k_effective_updater_i.h"

namespace bart {

namespace eigenvalue {

namespace k_effective {

class UpdaterViaFissionSourceI : public K_EffectiveUpdaterI {
 public:
  virtual ~UpdaterViaFissionSourceI() = default;
};

} // namespace k_effective

} // namespace eigenvalue

} // namespace bart

#endif // BART_SRC_EIGENVALUE_K_EFFECTIVE_UPDATER_VIA_FISSION_SOURCE_I_H_