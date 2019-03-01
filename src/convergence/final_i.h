#ifndef BART_SRC_CONVERGENCE_FINAL_I_H_
#define BART_SRC_CONVERGENCE_FINAL_I_H_

#include <optional>

#include "convergence/status.h"

namespace bart {

namespace convergence {

/*! \brief Interface for classes that check for final convergence of a
 * framework.
 *
 * Generally, it will check for the status of some system parameter (such as
 * flux or eigenvalue), or that a maximum number of iterations has been reached.
 *
 */

class FinalI {
 public:
  //! Typedef for value used for indexing iterations
  using IterationNumber = int;

  virtual ~FinalI() = default;

  /*! \brief Check for final convergence of the system.
   * In the general case, a side effect is that this is what will increment the
   * iteration count and updates the internal convergence::Status object.
   * \return a convergence::Status struct that contains the status of the
   * system convergence
   */
  virtual Status CheckFinalConvergence() = 0;

  /*! \brief Get status of system convergence
   *
   * \return a convergence::Status struct that contains the status of system
   * convergence.
   */
  virtual Status convergence_status() const = 0;

  /*! \brief Returns status of system completion.
   *
   * \return bool indicating status of system completion
   */
  virtual bool   convergence_is_complete() const = 0;

  /*! \brief Returns the set maximum iterations
   *
   * \return IterationNumber indicating maximum iterations.
   */
  virtual IterationNumber max_iterations() const = 0;

  /*! \brief Returns the current iteration.
   *
   * \return IterationNumber indicating current iteration.
   */
  virtual IterationNumber iteration() const = 0;

  /*! \brief Sets the maximum iteration.
   *
   * \param to_set value to set as maximum iterations.
   * \return reference to this object.
   */
  virtual FinalI& SetMaxIterations(IterationNumber to_set) = 0;

  /*! \brief Set the current iteration.
   *
   * \param to_set value to set as the current iteration.
   * \return reference to this object.
   */
  virtual FinalI& SetIteration(IterationNumber to_set) = 0;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_I_H_
