#ifndef BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_
#define BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_

#include "convergence/final.h"

namespace bart {

namespace convergence {

class FinalOrN : public Final {
 public:
  FinalOrN() = default;
  ~FinalOrN() = default;

  Status CheckFinalConvergence() override;

};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_