#ifndef BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_
#define BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_

#include "convergence/convergence_checker.hpp"

#include <deal.II/lac/vector.h>

//! Convergence checkers for problem vector moments
namespace bart::convergence::moments {

/*! \brief Checks for convergence between flux moments using the percentage
 * change in the L1 norm, compared to the current iteration (\f$i\f$):
 *
 * \f[
 *
 * \Delta_i = \frac{|\phi_i - \phi_{i-1}|_{1}}{|\phi_{i}|_{1}}
 *
 * \f]
 *
 * Convergence is achieved if \f$\Delta_i \leq \Delta_{\text{max}}\f$.
 * */

class ConvergenceCheckerL1Norm : public ConvergenceChecker<dealii::Vector<double>> {
 public:
  using Vector = dealii::Vector<double>;
  /*! \brief Default constructor, setting max delta to \f$10^{-6}\f$. */
  explicit ConvergenceCheckerL1Norm(const double max_delta = 1e-6);
  auto SetMaxDelta(const double& to_set) -> void override;
  auto IsConverged(const Vector& current_iteration, const Vector& previous_iteration) -> bool override;
};

} // namespace bart::convergence::moments

#endif // BART_SRC_CONVERGENCE_MOMENTS_SINGLE_MOMENT_L1_NORM_H_
