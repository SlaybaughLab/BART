#ifndef BART_SRC_CONVERGENCE_FINAL_H_
#define BART_SRC_CONVERGENCE_FINAL_H_

#include <optional>

#include "convergence/final_i.h"
#include "convergence/status.h"

namespace bart {

namespace convergence {

/*! \brief Implements getters and setters for final convergence interface.
 *
 */
template <typename CheckerType>
class Final : public FinalI<CheckerType> {
 public:
  using typename FinalI<CheckerType>::IterationNumber;

  virtual ~Final() = default;

  Status convergence_status() const override {
      return convergence_status_; };

  bool   convergence_is_complete() const override {
    return convergence_status_.is_complete; };

  IterationNumber max_iterations() const override {
      return convergence_status_.max_iterations; };

  IterationNumber iteration() const override {
      return convergence_status_.iteration_number; };

  Final<CheckerType>& SetMaxIterations(IterationNumber to_set) override;

  Final<CheckerType>& SetIteration(IterationNumber to_set) override;

 protected:
  Status convergence_status_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_H_
