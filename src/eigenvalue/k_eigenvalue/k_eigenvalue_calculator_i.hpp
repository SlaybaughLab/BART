#ifndef BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EIGENVALUE_CALCULATOR_I_HPP_
#define BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EIGENVALUE_CALCULATOR_I_HPP_

#include <optional>

#include "system/system.hpp"
#include "utility/has_description.h"

namespace bart::eigenvalue::k_eigenvalue {

/*! \brief Interface for classes that calculates the value of the _k_-eigenvalue.
 *
 * The _k_-eigenvalue problem is a version of the criticality problem in which the fission source is tuned by adjusting
 * the value of \f$\nu\f$ by a factor \f$k\f$ to achieve a time-independent result. The value of \f$k\f$ will then
 * identify if the fission source is too small, too large, or exactly correct to maintain criticality.
 *
 */
class K_EigenvalueCalculatorI : public utility::HasDescription {
 public:
  virtual ~K_EigenvalueCalculatorI() = default;

  /*! \brief Calculate the _k_-eigenvalue for the system.
   *
   * @param system system to calculate _k_ for.
   * @return double value of _k_.
   */
  virtual auto CalculateK_Eigenvalue(system::System& system) -> double = 0;

  /*! \brief Get the last value of _k_ calculated or std::nullopt if it has not been calculated yet. */
  virtual auto k_eigenvalue() const -> std::optional<double> = 0;
};

} // namespace bart::eigenvalue::k_eigenvalue

#endif // BART_SRC_EIGENVALUE_K_EIGENVALUE_K_EIGENVALUE_CALCULATOR_I_HPP_