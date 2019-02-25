#ifndef BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_
#define BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_

#include <memory>

#include "flux_checker_i.h"
#include "group_flux_checker_i.h"


namespace bart {

namespace convergence {

/*! \brief Checks each flux group sequentially for convergence.
 * Will not re-check converged groups until final convergence is believed to
 * have been reached. */

class GroupFluxCheckerSequential : public GroupFluxCheckerI {
 public:
  GroupFluxCheckerSequential() = default;
  GroupFluxCheckerSequential(std::unique_ptr<FluxCheckerI> &tester);
  ~GroupFluxCheckerSequential() = default;
  
  void ProvideChecker(std::unique_ptr<FluxCheckerI> &tester) {
    tester_ = std::move(tester); };
  
  bool isConverged(data::GroupFluxPointers &current,
                   data::GroupFluxPointers &last);
  bool isConverged() const override { return converged_; };
  int GetFailedGroup() const { return failed_group_; };

 private:
  /*! Flux convergence tester that will be used to check sequentially */
  std::unique_ptr<FluxCheckerI> tester_ = nullptr;
  /*! Stores the ID of the first group that failed the convergence check */
  int failed_group_ = 0;
  /*! Stores value of last isConverged call */
  bool converged_ = false;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_GROUP_FLUX_CHECKER_SEQUENTIAL_H_
