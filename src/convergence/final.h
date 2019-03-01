#ifndef BART_SRC_CONVERGENCE_FINAL_H_
#define BART_SRC_CONVERGENCE_FINAL_H_

#include <optional>

#include "convergence/final_i.h"
#include "convergence/status.h"

namespace bart {

namespace convergence {

class Final : public FinalI {
 public:
  Final() = default;
  virtual ~Final() = default;

  Status convergence_status() const override { return convergence_status_; };
  bool   convergence_is_complete() const override {
    return convergence_is_complete_; };
  IterationNumber max_iterations() const override { return max_iterations_; };
  IterationNumber iteration() const override { return iteration_; };
  void SetMaxIterations(IterationNumber to_set) override;
  void SetIteration(IterationNumber to_set) override;
 protected:
  Status convergence_status_;
  bool convergence_is_complete_ = false;
  IterationNumber max_iterations_ = 1;
  IterationNumber iteration_ = 0;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_H_
