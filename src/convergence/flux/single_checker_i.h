#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_

#include <optional>

#include "data/vector_parameters.h"

namespace bart {

namespace convergence {
 
namespace flux {

/*! \brief Checks for convergence between two provided fluxes.
 * Convergence is determined by calculating a delta between the two fluxes,
 * (generally using norms) and comparing them to a maximum allowed delta.
 */

class SingleCheckerI {
 public:
  virtual ~SingleCheckerI() = default;
  /* \brief Checks for convergence of two provided fluxes */
  virtual bool CheckIfConverged(data::FluxVector &, data::FluxVector &) = 0;
  /* \brief Returns status of convergence (from last call to CheckIfConverged */
  virtual bool is_converged() const = 0;
  /* \brief Set the threshold value for convergence check */
  virtual void SetMaxDelta(double to_set) = 0;
  /* \brief Get the threshold value for convergence check */
  virtual double max_delta() const = 0;
  /* \brief Get the delta value from the last convergence check. May be empty
   * if convergence has not be checked. */
  virtual std::optional<double> delta() const = 0;
};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_I_H_
