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
  int    max_iterations() const override { return max_iterations_; };
  void   SetMaxIterations(int to_set) override { max_iterations_ = to_set; };
 protected:
  Status convergence_status_;
  bool   convergence_is_complete_ = false;
  int    max_iterations_ = 1;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_H_
