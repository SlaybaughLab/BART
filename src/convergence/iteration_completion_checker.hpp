#ifndef BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_

#include <optional>

#include <deal.II/base/exceptions.h>

#include "convergence/iteration_completion_checker_i.hpp"
#include "convergence/status.hpp"

namespace bart::convergence {

/*! \brief Default implementation for the iterative convergence checker
 *
 */
template <typename CompareType>
class IterationCompletionChecker : public IterationCompletionCheckerI<CompareType> {
 public:
  using typename IterationCompletionCheckerI<CompareType>::IterationNumber;
  virtual ~IterationCompletionChecker() = default;

  [[nodiscard]] auto convergence_status() const -> Status override { return convergence_status_; };

  [[nodiscard]] auto convergence_is_complete() const -> bool override { return convergence_status_.is_complete; };

  [[nodiscard]] auto max_iterations() const -> IterationNumber override { return convergence_status_.max_iterations; };

  [[nodiscard]] auto iteration() const -> IterationNumber override { return convergence_status_.iteration_number; };

  auto SetMaxIterations(IterationNumber to_set) -> IterationCompletionChecker<CompareType>& override {
    AssertThrow(to_set > 0, dealii::ExcMessage("Max iterations must be > 0"))
    convergence_status_.max_iterations = to_set;
    return *this;
  }

  auto SetIteration(IterationNumber to_set) -> IterationCompletionChecker<CompareType>& override {
    AssertThrow(to_set >= 0, dealii::ExcMessage("Iteration must be >= 0"));
    convergence_status_.iteration_number = to_set;
    return *this;
  }

  auto Reset() -> void override {
    Status convergence_status;
    convergence_status.max_iterations = convergence_status_.max_iterations;
    convergence_status_ = convergence_status;
  };

 protected:
  Status convergence_status_;
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_
