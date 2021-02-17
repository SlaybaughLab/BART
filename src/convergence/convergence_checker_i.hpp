#ifndef BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_I_HPP_
#define BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_I_HPP_

#include <optional>

//! Tools for checking the convergence of variables
namespace bart::convergence {

/*! \brief Checks for convergence between two provided values.
 * Convergence is determined by calculating a delta between the two values,
 * and comparing them to a maximum allowed delta.
 *
 * \tparam CompareT type of value to compare.
 * \tparam DeltaT type of delta that will be used to determine comparison, default double.
 */

template <typename CompareT, typename DeltaT>
class ConvergenceCheckerI {
 public:
  virtual ~ConvergenceCheckerI() = default;
  /*! \brief Checks for convergence of two provided values */
  virtual auto CheckIfConverged(const CompareT& current_iteration, const CompareT& previous_iteration) -> bool = 0;
  /*! \brief Returns status of convergence (from last call to CheckIfConverged */
  virtual auto is_converged() const -> bool = 0;
  /*! \brief Set the threshold value for convergence check */
  virtual auto SetMaxDelta(const DeltaT& to_set) -> void = 0;
  /*! \brief Get the threshold value for convergence check */
  virtual auto max_delta() const -> DeltaT = 0;
  /*! \brief Get the delta value from the last convergence check. May be empty
   * if convergence has not be checked. */
  virtual auto delta() const -> std::optional<DeltaT> = 0;
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_I_HPP_