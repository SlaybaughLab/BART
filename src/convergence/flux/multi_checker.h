#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_

#include <memory>
#include <optional>

#include "../../data/vector_parameters.h"
#include "multi_checker_i.h"
#include "single_checker_i.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Checks that all fluxes have converged */

class MultiChecker : public MultiCheckerI{
 public:
  ~MultiChecker() = default;
  void ProvideChecker(std::unique_ptr<SingleCheckerI> &checker) {
    checker_ = std::move(checker); };
  bool is_converged() const override { return is_converged_; }
  std::optional<int> GetFailedIndex() const override;
 protected:
  std::unique_ptr<SingleCheckerI> checker_ = nullptr;
  bool is_converged_ = false;
  int failed_index_ = 0;
};

std::optional<int> MultiChecker::GetFailedIndex() const {
  if (!is_converged_)
    return failed_index_;
  else
    return {};
}

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_
