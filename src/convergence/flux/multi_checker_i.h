#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_

#include <memory>
#include <optional>

#include "data/vector_parameters.h"
#include "convergence/flux/single_checker_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks that all fluxes have converged */

class MultiCheckerI {
 public:
  virtual ~MultiCheckerI() = default;
  /*! \brief Check for convergence of all fluxes provided */
  virtual bool CheckIfConverged(data::MultiFluxPtrs &current_iteration,
                                data::MultiFluxPtrs &previous_iteration) = 0;
  
  /* \brief Returns status of convergence (from last call to CheckIfConverged) */
  virtual bool is_converged() const = 0;

  /* \brief Returns the index of the flux that failed convergence check */
  virtual std::optional<int> failed_index() const = 0;

  /* \brief Returns the delta of the flux that failed convergence check */
  virtual std::optional<double> failed_delta() const = 0;
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_I_H_