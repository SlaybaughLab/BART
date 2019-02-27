#ifndef BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_
#define BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_

#include <optional>

#include "single_checker_i.h"
#include "../../data/vector_parameters.h"

namespace bart {

namespace convergence {
 
namespace flux {

/*! \brief Checks for convergence between two provided fluxes. */

class SingleChecker : public SingleCheckerI {
 public:
  SingleChecker() = default;
  explicit SingleChecker(double max_delta) : max_delta_(max_delta) {}
  virtual ~SingleChecker() = default;
  bool is_converged() const override { return is_converged_; };
  void SetMaxDelta(double to_set) override { max_delta_ = to_set; };
  double max_delta() const override { return max_delta_; }
  std::optional<double> delta() const { return delta_; };
 protected:
  /*! Delta between fluxes from last convergence check */
  std::optional<double> delta_ = std::nullopt;
  
  /*! Did last convergence check result in convergence */
  bool is_converged_ = false;
  
  /*! Maximum delta for convergence */
  double max_delta_ = 0;

};

} // namespace flux
  
} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_SINGLE_CHECKER_H_
