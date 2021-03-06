#ifndef BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_HPP_

#include <optional>

#include <deal.II/base/exceptions.h>

#include "convergence/convergence_checker_i.hpp"

namespace bart::convergence {

/* This concept ensures that the type we are using for delta can be default constructed for the initial value of
 * max_delta_, and copy assignable to ensure SetMaxDelta doesn't cause issues */
template <typename T>
concept ConstructableAndCopyable = std::is_default_constructible_v<T> && std::is_copy_assignable_v<T>;

/*! \brief Default implementation of ConvergenceCheckerI.
 *
 * Checks for the convergence of two values, based on a maximum delta. No actual check is defined in this class,
 * as it doesn't even know how to calculate delta (or in some cases what its time is).
 *
 * @tparam CompareT type to be compared.
 * @tparam DeltaT type for delta (default double), but it can be anything that is copy assignable (to be stored
 * internally) and default constructible (for initialization).
 */
template <typename CompareT, ConstructableAndCopyable DeltaT = double>
class ConvergenceChecker : public ConvergenceCheckerI<CompareT, DeltaT> {
 public:
  virtual ~ConvergenceChecker() = default;

  [[nodiscard]] auto is_converged() const -> bool override { return is_converged_; };
  auto SetMaxDelta(const DeltaT& to_set) -> void override { max_delta_ = to_set; };
  [[nodiscard]] auto max_delta() const -> DeltaT override { return max_delta_; }
  [[nodiscard]] auto delta() const -> std::optional<DeltaT> override { return delta_; };
  [[nodiscard]] auto failed_index() const -> std::optional<int> override { return failed_index_; };
 protected:
  /*! Delta between values from last convergence check */
  std::optional<DeltaT> delta_{ std::nullopt };
  /*! Failed index from last convergence check */
  std::optional<int> failed_index_{ std::nullopt };
  /*! Did last convergence check result in convergence */
  bool is_converged_{ false };
  /*! Maximum delta for convergence */
  DeltaT max_delta_{ DeltaT() };
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_CONVERGENCE_CHECKER_HPP_