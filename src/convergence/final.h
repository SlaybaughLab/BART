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
template <typename CompareType>
class Final : public FinalI<CompareType> {
 public:
  using typename FinalI<CompareType>::IterationNumber;

  virtual ~Final() = default;

  Status convergence_status() const override {
      return convergence_status_; };

  bool convergence_is_complete() const override {
    return convergence_status_.is_complete; };

  IterationNumber max_iterations() const override {
      return convergence_status_.max_iterations; };

  IterationNumber iteration() const override {
      return convergence_status_.iteration_number; };

  Final<CompareType>& SetMaxIterations(IterationNumber to_set) override;

  Final<CompareType>& SetIteration(IterationNumber to_set) override;

  void Reset() override {
    Status convergence_status;
    convergence_status.max_iterations = convergence_status_.max_iterations;
    convergence_status_ = convergence_status;
  };

 protected:
  Status convergence_status_;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FINAL_H_
