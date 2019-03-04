#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_

#include <memory>
#include <optional>

#include "data/vector_parameters.h"
#include "convergence/flux/single_checker_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks that all fluxes have converged.
 *
 * Fluxes are passed via data::MultiFluxPtrs, for the current and previous
 * iteration.
 * */
class MultiCheckerI {
 public:
  virtual ~MultiCheckerI() = default;

  /*! \brief Identifies if all fluxes provided have converged.
   *
   * \param current_iteration fluxes for the current iteration
   * \param previous_iteration fluxes for the previous iteration
   * \return bool indicating status of convergence
   */
  virtual bool CheckIfConverged(data::ScalarFluxPtrs &current_iteration,
                                data::ScalarFluxPtrs &previous_iteration) = 0;
  
  /*! \brief Returns status of previous call to CheckIfConverged
   *
   * \return bool indicating if convergence is reached.
   */
  virtual bool is_converged() const = 0;

  /*! \brief Identifies the index that caused the convergence check to fail.
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

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_
