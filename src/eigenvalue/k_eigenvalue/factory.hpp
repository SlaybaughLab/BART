#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_FACTORY_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_FACTORY_HPP_

namespace bart::eigenvalue::k_eigenvalue {

enum class K_EffectiveUpdaterName {
  kCalculatorViaFissionSource = 0,
  kCalculatorViaRayleighQuotient = 1
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif //BART_SRC_EIGENVALUE_K_EIGENVALUE_FACTORY_HPP_
