#ifndef BART_SRC_CONVERGENCE_SINGLE_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_SINGLE_CHECKER_HPP_

#include <optional>

#include <deal.II/base/exceptions.h>

#include "convergence/single_checker_i.hpp"

namespace bart::convergence {

/* This concept ensures that the type we are using for delta can be default constructed for the initial value of
 * max_delta_, and copy assignable to ensure SetMaxDelta doesn't cause issues */
template <typename T>
concept ConstructableAndCopyable = std::is_default_constructible_v<T> && std::is_copy_assignable_v<T>;

template <typename CompareT, ConstructableAndCopyable DeltaT = double>
class SingleChecker : public SingleCheckerI<CompareT, DeltaT> {
 public:
  virtual ~SingleChecker() = default;

  [[nodiscard]] auto is_converged() const -> bool override { return is_converged_; };
  auto SetMaxDelta(const DeltaT& to_set) -> void override { max_delta_ = to_set; };
  [[nodiscard]] auto max_delta() const -> DeltaT override { return max_delta_; }
  [[nodiscard]] auto delta() const -> std::optional<DeltaT> override { return delta_; };
 protected:
  /*! Delta between moments from last convergence check */
  std::optional<DeltaT> delta_{ std::nullopt };
  /*! Did last convergence check result in convergence */
  bool is_converged_{ false };
  /*! Maximum delta for convergence */
  DeltaT max_delta_{ DeltaT() };
};

} // namespace bart::convergence

#endif // BART_SRC_CONVERGENCE_SINGLE_CHECKER_HPP_