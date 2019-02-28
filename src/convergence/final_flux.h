#ifndef BART_SRC_CONVERGENCE_FINAL_FLUX_H_
#define BART_SRC_CONVERGENCE_FINAL_FLUX_H_

#include <memory>

#include "convergence/final.h"
#include "convergence/status.h"
#include "convergence/flux/multi_checker_i.h"
#include "utility/uncopyable.h"

/*! \brief Checks for final convergence of flux, or max iterations reached. */

namespace bart {

namespace convergence {

class FinalFlux : public Final, private utility::Uncopyable {
 public:
  FinalFlux() = default;
  ~FinalFlux() = default;

  void ProvideMultiChecker(std::unique_ptr<flux::MultiCheckerI> &checker);
  Status CheckFinalConvergence() override;
 private:
  std::unique_ptr<flux::MultiCheckerI> checker_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_FLUX_H_
