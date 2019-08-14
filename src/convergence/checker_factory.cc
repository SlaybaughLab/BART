#include "convergence/checker_factory.h"

#include "convergence/moments/single_moment_checker_l1_norm.h"
#include "convergence/parameters/single_parameter_checker.h"
#include "convergence/final_checker_or_n.h"
#include "convergence/final_i.h"

namespace bart {

namespace convergence {

std::unique_ptr<FinalI<double>> CheckerFactory::MakeParameterChecker(
    const double max_delta,
    const int max_iterations,
    ParameterCheckerImpl parameter_checker_type,
    FinalCheckerImpl final_checker_type) {
  std::unique_ptr<FinalI<double>> return_ptr = nullptr;

  if (final_checker_type == FinalCheckerImpl::kConvergenceOrN) {
    if (parameter_checker_type == ParameterCheckerImpl::kDefault) {

      auto checker_ptr = std::make_unique<parameters::SingleParameterChecker>(
          max_delta);

      using ReturnType = FinalCheckerOrN<double,
                                         parameters::SingleParameterChecker>;

      return_ptr = std::move(
          std::make_unique<ReturnType>(std::move(checker_ptr)));
      return_ptr->SetMaxIterations(max_iterations);
    }
  }

  return return_ptr;
}

std::unique_ptr<FinalI<system::moments::MomentVector>>
CheckerFactory::MakeSingleMomentChecker(
    const double max_delta,
    const int max_iterations,
    SingleMomentCheckerImpl single_moment_checker_type,
    FinalCheckerImpl final_checker_type) {

  std::unique_ptr<FinalI<system::moments::MomentVector>> return_ptr = nullptr;

  if (final_checker_type == FinalCheckerImpl::kConvergenceOrN) {
    if (single_moment_checker_type == SingleMomentCheckerImpl::kL1Norm) {
      auto checker_ptr = std::make_unique<moments::SingleMomentCheckerL1Norm>(
          max_delta);

      using ReturnType = FinalCheckerOrN<system::moments::MomentVector,
                                         moments::SingleMomentCheckerI>;

      return_ptr = std::move(
          std::make_unique<ReturnType>(std::move(checker_ptr)));
      return_ptr->SetMaxIterations(max_iterations);
    }
  }

  return return_ptr;
}
} // namespace convergence

} //namespace bart

