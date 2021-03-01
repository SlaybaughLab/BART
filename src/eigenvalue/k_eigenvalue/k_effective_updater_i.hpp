#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EFFECTIVE_UPDATER_I_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EFFECTIVE_UPDATER_I_HPP_

#include <optional>

#include "system/system.hpp"
#include "utility/has_description.h"

namespace bart::eigenvalue::k_eigenvalue {

class K_EffectiveUpdaterI : public utility::HasDescription {
 public:
  virtual ~K_EffectiveUpdaterI() = default;

  virtual auto k_effective() const -> std::optional<double> = 0;
  virtual auto CalculateK_Effective(system::System& system) -> double = 0;
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EFFECTIVE_UPDATER_I_HPP_