#ifndef BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_
#define BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_

#include <memory>

#include "convergence/final.h"
#include "convergence/status.h"

namespace bart {

namespace convergence {

template <typename CheckerType>
class FinalCheckerOrN : public Final{
 public:
  FinalCheckerOrN(std::unique_ptr<CheckerType> checker_ptr)
      : checker_ptr_(std::move(checker_ptr)) {}
  ~FinalCheckerOrN() = default;

  Status CheckFinalConvergence() override;

 protected:
  std::unique_ptr<CheckerType> checker_ptr_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_