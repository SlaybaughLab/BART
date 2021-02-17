#ifndef BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_I_HPP_
#define BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_I_HPP_

#include <optional>

#include "convergence/status.hpp"

namespace bart::convergence {

/*! \brief Interface for classes that determines the status of final
 * convergence.
 *
 * Generally, it will check for the status of some system parameter (such as
 * flux or eigenvalue), or that a maximum number of iterations has been reached.
 *
 * \tparam CompareType the types of the objects or values that will be compared
 * to determine final convergence. (i.e. fluxes, integers, etc).
 *
 */
template <typename CompareType>
class IterationCompletionCheckerI {
 public:
  //! Typedef for value used for indexing iterations
  using IterationNumber = int;

  virtual ~IterationCompletionCheckerI() = default;

  /*! \brief Check for final convergence of the system.
   * In the general case, a side effect is that this is what will increment the
   * iteration count and updates the internal convergence::Status object.
   * \return a convergence::Status struct that contains the status of the
   * system convergence
   */
  virtual auto ConvergenceStatus(CompareType& current_iteration, CompareType& previous_iteration) -> Status = 0;

  /*! \brief Get status of system convergence
   *
   * \return a convergence::Status struct that contains the status of system
   * convergence.
   */
  virtual auto convergence_status() const -> Status = 0;

  /*! \brief Returns status of system completion.
   *
   * \return bool indicating status of system completion
   */
  virtual auto convergence_is_complete() const -> bool = 0;

  /*! \brief Returns the set maximum iterations
   *
   * \return IterationNumber indicating maximum iterations.
   */
  virtual auto max_iterations() const -> IterationNumber = 0;

  /*! \brief Returns the current iteration.
   *
   * \return IterationNumber indicating current iteration.
   */
  virtual auto iteration() const -> IterationNumber = 0;

  /*! \brief Sets the maximum iteration.
   *
   * \param to_set value to set as maximum iterations.
   * \return reference to this object.
   */
  virtual auto SetMaxIterations(IterationNumber to_set) -> IterationCompletionCheckerI& = 0;

  /*! \brief Set the current iteration.
   *
   * \param to_set value to set as the current iteration.
   * \return reference to this object.
   */
  virtual auto SetIteration(IterationNumber to_set) -> IterationCompletionCheckerI& = 0;

  /*! \brief Resets the iterative convergence checker */
  virtual auto Reset() -> void = 0;
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_I_HPP_
