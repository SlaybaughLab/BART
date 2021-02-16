#ifndef BART_CONVERGENCE_SINGLE_CHECKER_I_H_
#define BART_CONVERGENCE_SINGLE_CHECKER_I_H_

#include <optional>

namespace bart {

namespace convergence {

/*! \brief Checks for convergence between two provided values.
 * Convergence is determined by calculating a delta between the two values,
 * and comparing them to a maximum allowed delta.
 *
 * \tparam CompareType type of value to compare.
 *
 */

template <typename CompareType>
class SingleCheckerI {
 public:
  virtual ~SingleCheckerI() = default;
  /*! \brief Checks for convergence of two provided values */
  virtual bool CheckIfConverged(const CompareType& current_iteration,
                                const CompareType& previous_iteration) = 0;
  /*! \brief Returns status of convergence (from last call to CheckIfConverged */
  virtual bool is_converged() const = 0;
  /*! \brief Set the threshold value for convergence check */
  virtual void SetMaxDelta(const double to_set) = 0;
  /*! \brief Get the threshold value for convergence check */
  virtual double max_delta() const = 0;
  /*! \brief Get the delta value from the last convergence check. May be empty
   * if convergence has not be checked. */
  virtual std::optional<double> delta() const = 0;
};



} // namespace convergence

} // namespace bart

#endif // BART_CONVERGENCE_SINGLE_CHECKER_I_H_