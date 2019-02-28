#ifndef BART_SRC_CONVERGENCE_FINAL_FLUX_H_
#define BART_SRC_CONVERGENCE_FINAL_FLUX_H_

#include "convergence/final.h"
#include "utility/uncopyable.h"

/*! \brief Checks for final convergence of flux, or max iterations reached. */

namespace bart {

namespace convergence {

class FinalFlux : public Final, private Uncopyable {
 public:
  FinalFlux() = default;
  ~FinalFlux() = default;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_FLUX_H_
