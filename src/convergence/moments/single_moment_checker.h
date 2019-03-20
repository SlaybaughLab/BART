#ifndef BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_H_
#define BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_H_

#include <optional>

#include "convergence/moments/single_moment_checker_i.h"
#include "data/moment_types.h"

namespace bart {

namespace convergence {
 
namespace moments {

/*! \brief Checks for convergence between two provided moments. */

class SingleMomentChecker : public SingleMomentCheckerI {
 public:
  SingleMomentChecker() = default;
  explicit SingleMomentChecker(const double max_delta)
      : max_delta_(max_delta) {}

  virtual ~SingleMomentChecker() = default;

  bool is_converged() const override { return is_converged_; };
  void SetMaxDelta(const double to_set) override { max_delta_ = to_set; };
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

} // namespace moments
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_H_
