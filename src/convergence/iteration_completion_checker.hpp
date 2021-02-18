#ifndef BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_

#include <memory>
#include <optional>

#include <deal.II/base/exceptions.h>

#include "convergence/convergence_checker_i.hpp"
#include "convergence/iteration_completion_checker_i.hpp"
#include "convergence/status.hpp"

namespace bart::convergence {

/*! \brief Default implementation for the iterative convergence checker
 *
 * This class will use a provided convergence checker to check for convergence and update a convergence::Status
 * object appropriately. Convergence is considered complete if the values passed have converged or max iterations
 * have been reached.
 *
 */
template <typename CompareType>
class IterationCompletionChecker : public IterationCompletionCheckerI<CompareType> {
 public:
  using ConvergenceChecker = typename convergence::ConvergenceCheckerI<CompareType>;
  using typename IterationCompletionCheckerI<CompareType>::IterationNumber;
  [[deprecated]] IterationCompletionChecker() = default;
  IterationCompletionChecker(std::unique_ptr<ConvergenceChecker>);
  virtual ~IterationCompletionChecker() = default;

  [[nodiscard]] auto ConvergenceStatus(CompareType &current_iteration, CompareType &previous_iteration) -> Status override;
  [[nodiscard]] auto convergence_status() const -> Status override { return convergence_status_; };
  [[nodiscard]] auto convergence_is_complete() const -> bool override { return convergence_status_.is_complete; };
  [[nodiscard]] auto max_iterations() const -> IterationNumber override { return convergence_status_.max_iterations; };
  [[nodiscard]] auto iteration() const -> IterationNumber override { return convergence_status_.iteration_number; };

  auto SetMaxIterations(IterationNumber to_set) -> IterationCompletionChecker<CompareType>& override;
  auto SetIteration(IterationNumber to_set) -> IterationCompletionChecker<CompareType>& override;
  auto Reset() -> void override;

  auto convergence_checker_ptr() -> ConvergenceChecker* { return convergence_checker_ptr_.get(); }
 protected:
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_{ nullptr };
  Status convergence_status_;
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_ITERATION_COMPLETION_CHECKER_HPP_
