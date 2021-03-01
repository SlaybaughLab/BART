#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_I_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_I_HPP_

#include "eigenvalue/k_eigenvalue/k_eigenvalue_calculator_i.hpp"

namespace bart::eigenvalue::k_eigenvalue {

class CalculatorViaFissionSourceI : public K_EigenvalueCalculatorI {
 public:
  virtual ~CalculatorViaFissionSourceI() = default;

  /*! \brief Return the current problem fission source. */
  virtual auto current_fission_source() const -> std::optional<double> = 0;
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_CALCULATOR_VIA_FISSION_SOURCE_I_HPP_