#ifndef BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_I_H_

#include <memory>
#include <optional>

#include "convergence/moments/single_moment_checker_i.h"
#include "data/moment_types.h"

namespace bart {

namespace convergence {

namespace moments {

/*! \brief Checks that all fluxes have converged.
 *
 * Fluxes are passed via a data::MomentsMap. Which moments are checked is up to
 * the implementation of the checker.
 * */
class MultiCheckerI {
 public:
  virtual ~MultiCheckerI() = default;

  /*! \brief Identifies if all moments provided have converged.
   *
   * \param current_iteration all moments for current iteration.
   * \param previous_iteration all moments for previous iteration
   * \return bool indicating if convergence has been reached.
   */
  virtual bool CheckIfConverged(const data::MomentsMap &current_iteration,
                                const data::MomentsMap &previous_iteration) = 0;
  
  /*! \brief Returns status of previous call to CheckIfConverged
   *
   * \return bool indicating if convergence is reached.
   */
  virtual bool is_converged() const = 0;

  /*! \brief Identifies the index that caused the convergence check to fail.
   *
   * Specific implementation determines the value that will be stored.
   *
   * \return failed index via std::optional<int>, will be empty if converged or
   * other error.
   */
  virtual std::optional<int> failed_index() const = 0;

  /*! \brief Returns the delta that resulted in failed convergence.
   *
   * \return delta via std::optional<double>, will be empty if no delta.
   */
  virtual std::optional<double> delta() const = 0;
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_I_H_
