#ifndef BART_SRC_CONVERGENCE_CHECKER_FACTORY_H_
#define BART_SRC_CONVERGENCE_CHECKER_FACTORY_H_

#include "system/moments/spherical_harmonic_types.h"

#include <memory>

namespace bart {

namespace convergence {
template <typename CompareType> class FinalI;

enum class SingleMomentCheckerImpl {
  kL1Norm = 0,
};

enum class ParameterCheckerImpl {
  kDefault = 0,
};

enum class FinalCheckerImpl {
  kConvergenceOrN = 0,
};

class CheckerFactory {
 public:
  static std::unique_ptr<FinalI<double>> MakeParameterChecker(
      const double max_delta = 1e-6,
      const int max_iterations = 100,
      ParameterCheckerImpl parameter_checker_type = ParameterCheckerImpl::kDefault,
      FinalCheckerImpl final_checker_type = FinalCheckerImpl::kConvergenceOrN);

  static std::unique_ptr<FinalI<system::moments::MomentVector>>
  MakeSingleMomentChecker(
      const double max_delta = 1e-6,
      const int max_iterations = 100,
      SingleMomentCheckerImpl single_moment_checker_type = SingleMomentCheckerImpl::kL1Norm,
      FinalCheckerImpl final_checker_type = FinalCheckerImpl::kConvergenceOrN);
};

} // namespace convergence

} //namespace bart

#endif //BART_SRC_CONVERGENCE_CHECKER_FACTORY_H_
