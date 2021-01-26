#ifndef BART_SRC_EIGENVALUE_K_EFFECTIVE_FACTORY_HPP_
#define BART_SRC_EIGENVALUE_K_EFFECTIVE_FACTORY_HPP_

namespace bart::eigenvalue::k_effective {

enum class K_EffectiveUpdaterName {
  kUpdaterViaFissionSource = 0,
  kUpdaterViaRayleighQuotient = 1
};

} // namespace bart::eigenvalue::k_effective

#endif //BART_SRC_EIGENVALUE_K_EFFECTIVE_FACTORY_HPP_
