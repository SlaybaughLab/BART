#ifndef BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_
#define BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_

#include "convergence/convergence_checker.hpp"
#include "system/moments/spherical_harmonic_types.h"

namespace bart {

namespace convergence {

namespace moments {

/*! \brief Checks for convergence between two provided moments.
 * Convergence is determined by calculating a delta between the two moments,
 * (generally using norms) and comparing them to a maximum allowed delta.
 */

class SingleMomentCheckerI : public ConvergenceChecker<system::moments::MomentVector> {
 public:
  virtual ~SingleMomentCheckerI() = default;

 protected:
  using ConvergenceChecker<system::moments::MomentVector>::max_delta_;
  using ConvergenceChecker<system::moments::MomentVector>::delta_;
  using ConvergenceChecker<system::moments::MomentVector>::is_converged_;
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_CHECKER_I_H_