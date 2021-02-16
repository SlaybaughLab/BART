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

  bool is_converged() const override { return is_converged_; };
  void SetMaxDelta(const DeltaT& to_set) override { max_delta_ = to_set; };
  DeltaT max_delta() const override { return max_delta_; }
  std::optional<DeltaT> delta() const { return delta_; };
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