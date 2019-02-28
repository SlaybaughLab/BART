#ifndef BART_SRC_CONVERGENCE_FINAL_FLUX_H_
#define BART_SRC_CONVERGENCE_FINAL_FLUX_H_

#include "convergence/final.h"
#include "convergence/status.h"
#include "utility/uncopyable.h"

/*! \brief Checks for final convergence of flux, or max iterations reached. */

namespace bart {

namespace convergence {

class FinalFlux : public Final, private utility::Uncopyable {
 public:
  FinalFlux() = default;
  ~FinalFlux() = default;

  Status CheckFinalConvergence() override { return convergence_status_; };
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_FLUX_H_