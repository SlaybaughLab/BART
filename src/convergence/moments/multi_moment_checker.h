#ifndef BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_H_
#define BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_H

#include "system/moments/moment_types.h"
#include "convergence/moments/multi_moment_checker_i.h"
#include "utility/uncopyable.h"

namespace bart {

namespace convergence {

namespace moments {

/*! \brief Implementation of most getting and setting base class funtions for
 * convergence::moments::MultiMomentCheckerI objects.
 */

class MultiMomentChecker : public MultiMomentCheckerI, private utility::Uncopyable {
 public:
  /*! \brief Constructor, taking a SingleChecker that will be used to compare
   * moments for convergence.
   *
   * \param checker to use.
   */
  explicit MultiMomentChecker(std::unique_ptr<SingleMomentCheckerI> checker)
      : checker_(std::move(checker)) {};

  ~MultiMomentChecker() = default;
  bool is_converged() const override { return is_converged_; }
  std::optional<int> failed_index() const override { return failed_index_; };
  std::optional<double> delta() const override { return delta_; };
 protected:
  std::unique_ptr<SingleMomentCheckerI> checker_;
  bool is_converged_ = false;
  std::optional<int> failed_index_ = std::nullopt;
  std::optional<double> delta_ = std::nullopt;
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_H_
