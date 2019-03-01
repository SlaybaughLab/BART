#ifndef BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_
#define BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_

#include <memory>

#include "convergence/final.h"
#include "convergence/status.h"
#include "convergence/flux/multi_checker_i.h"
#include "data/system_fluxes.h"
#include "utility/uncopyable.h"

namespace bart {

namespace convergence {

/*! \brief Checks for final convergence of flux, or max iterations reached.
 *
 * Requires access to a
 * \ref flux::MultiCheckerI, which it will use to
 * compare the current and previous flux iterations contained in a
 * \ref bart::data::SystemFluxes object. Will also return a state of converged if
 * iterations are reached.
*
 * \param checker a pointer to a \ref bart::convergence::flux::MultiCheckerI
 * object, this class takes ownership of the object.
 *
 * \param fluxes a shared pointer to the system fluxes.
 *
 */
class FinalFluxOrN : public Final, private utility::Uncopyable {
 public:
  FinalFluxOrN(std::unique_ptr<flux::MultiCheckerI> checker,
            std::shared_ptr<data::SystemFluxes> fluxes);
  ~FinalFluxOrN() = default;

  Status CheckFinalConvergence() override;
 private:
  std::unique_ptr<flux::MultiCheckerI> checker_;
  std::shared_ptr<data::SystemFluxes>  fluxes_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_FLUX_OR_N_H_
