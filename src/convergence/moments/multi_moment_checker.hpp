#ifndef BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_HPP_
#define BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_HPP_

#include "convergence/convergence_checker.hpp"
#include "system/moments/spherical_harmonic_types.h"
#include "utility/uncopyable.h"

namespace bart::convergence::moments {


/*! \brief Implementation of most getting and setting base class funtions for
 * convergence::moments::MultiMomentCheckerI objects.
 */
class MultiMomentChecker : public ConvergenceChecker<system::moments::MomentsMap>, private utility::Uncopyable {
 public:
  using SingleMomentChecker = convergence::ConvergenceCheckerI<system::moments::MomentVector>;
  /*! \brief Constructor, taking a SingleChecker that will be used to compare
   * moments for convergence.
   *
   * \param checker to use.
   */
  explicit MultiMomentChecker(std::unique_ptr<SingleMomentChecker> checker) : checker_(std::move(checker)) {};

  virtual ~MultiMomentChecker() = default;
  auto is_converged() const -> bool override { return is_converged_; }
  auto failed_index() const -> std::optional<int> override { return failed_index_; };
  auto delta() const -> std::optional<double> override { return delta_; };
 protected:
  std::unique_ptr<SingleMomentChecker> checker_{ nullptr };
  bool is_converged_{ false };
  std::optional<int> failed_index_{ std::nullopt };
  std::optional<double> delta_{ std::nullopt };
};

} // namespace bart::convergence::moments

#endif // BART_SRC_CONVERGENCE_CONVERGENCE_MOMENTS_MULTI_MOMENT_CHECKER_HPP_
