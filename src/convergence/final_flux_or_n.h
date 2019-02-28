#ifndef BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_
#define BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_

#include <memory>

#include "convergence/final.h"
#include "convergence/status.h"
#include "convergence/flux/multi_checker_i.h"
#include "data/forward_declarations.h"
#include "utility/uncopyable.h"

namespace bart {

namespace convergence {

/*! \brief Checks for final convergence of flux, or max iterations reached.
 *
 * Requires access to a convergence::flux::MultiCheckerI, which it will use to
 * compare the current and previous flux iterations contained in an object of 
 * type data::SystemFluxes. Will also return a state of converged if max
 * iterations are reached.
 *
 */
class FinalFluxOrN : public Final, private utility::Uncopyable {
 public:
  FinalFluxOrN(std::unique_ptr<flux::MultiCheckerI> checker,
            std::shared_ptr<data::SystemFluxes> fluxes)
      : checker_(std::move(checker)),
        fluxes_(fluxes) {};
  ~FinalFluxOrN() = default;

  Status CheckFinalConvergence() override;
 private:
  std::unique_ptr<flux::MultiCheckerI> checker_;
  std::shared_ptr<data::SystemFluxes>  fluxes_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_
