#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_

#include <memory>
#include <optional>

#include "convergence/flux/multi_checker_i.h"
#include "convergence/flux/single_checker_i.h"
#include "data/vector_parameters.h"
#include "utility/uncopyable.h"

namespace bart {

namespace convergence {

namespace flux {

/*! \brief Implementation of most getting and setting base class funtions for
 * convergence::flux::MultiCheckerI objects.
 */

class MultiChecker : public MultiCheckerI, private utility::Uncopyable {
 public:
  explicit MultiChecker(std::unique_ptr<SingleCheckerI> checker)
      : checker_(std::move(checker)) {};
  ~MultiChecker() = default;
  bool is_converged() const override { return is_converged_; }
  std::optional<int> failed_index() const override { return failed_index_; };
  std::optional<double> delta() const override { return delta_; };
 protected:
  std::unique_ptr<SingleCheckerI> checker_;
  bool is_converged_ = false;
  std::optional<int> failed_index_ = std::nullopt;
  std::optional<double> delta_ = std::nullopt;
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_MULTI_CHECKER_H_
