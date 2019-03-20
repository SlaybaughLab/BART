#ifndef BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_

#include "data/moment_types.h"

namespace bart {

namespace data {

namespace moments {

/*! \brief Checks for convergence between two provided moments.
 * Convergence is determined by calculating a delta between the two moments,
 * (generally using norms) and comparing them to a maximum allowed delta.
 */

class SingleMomentCheckerI {
 public:
  virtual ~SingleMomentCheckerI() = default;
  /* \brief Checks for convergence of two provided fluxes */
  virtual bool CheckIfConverged(const data::MomentVector &,
                                const data::MomentVector&) = 0;
  /* \brief Returns status of convergence (from last call to CheckIfConverged */
  virtual bool is_converged() const = 0;
  /* \brief Set the threshold value for convergence check */
  virtual void SetMaxDelta(const double to_set) = 0;
  /* \brief Get the threshold value for convergence check */
  virtual double max_delta() const = 0;
  /* \brief Get the delta value from the last convergence check. May be empty
   * if convergence has not be checked. */
  virtual std::optional<double> delta() const = 0;
};

} // namespace moments

} // namespace data

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_