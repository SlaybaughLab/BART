#ifndef BART_SRC_CONVERGENCE_SINGLE_CHECKER_H_
#define BART_SRC_CONVERGENCE_SINGLE_CHECKER_H_

#include <optional>

#include "deal.II/base/exceptions.h"

#include "convergence/single_checker_i.h"

namespace bart {

namespace convergence {

template <typename CompareType>
class SingleChecker : public SingleCheckerI<CompareType> {
 public:
  virtual ~SingleChecker() = default;

  bool is_converged() const override { return is_converged_; };
  void SetMaxDelta(const double to_set) override {
    AssertThrow(to_set >= 0,
        dealii::ExcMessage("Max delta value must be >= 0"));
    max_delta_ = to_set; };
  double max_delta() const override { return max_delta_; }
  std::optional<double> delta() const { return delta_; };

 protected:
  /*! Delta between moments from last convergence check */
  std::optional<double> delta_ = std::nullopt;

  /*! Did last convergence check result in convergence */
  bool is_converged_ = false;

  /*! Maximum delta for convergence */
  double max_delta_ = 0;

};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_SINGLE_CHECKER_H_