#ifndef BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_
#define BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_

#include <memory>

#include "convergence/final.h"
#include "convergence/status.h"
#include "convergence/reporter/mpi_i.h"

namespace bart {

namespace convergence {

/*! \brief Check convergence via a checker class, or N iterations.
 *
 * \tparam CompareType the types of objects that will be compared to determine
 * convergence.
 * \tparam CheckerType type of checker used to determine convergence.
 */

template <typename CompareType, typename CheckerType>
class FinalCheckerOrN : public Final<CompareType>{
 public:
  /*! \brief Constructor.
   *
   * \param checker_ptr pointer to the convergence checker that this class
   * will take ownership of.
   */
  explicit FinalCheckerOrN(std::unique_ptr<CheckerType> checker_ptr)
      : FinalCheckerOrN(std::move(checker_ptr), nullptr) {}
  FinalCheckerOrN(std::unique_ptr<CheckerType> checker_ptr,
      std::shared_ptr<reporter::MpiI> reporter_ptr)
      : checker_ptr_(std::move(checker_ptr)),
        reporter_ptr_(reporter_ptr) {}
  ~FinalCheckerOrN() = default;

  Status CheckFinalConvergence(CompareType& current_iteration,
                               CompareType& previous_iteration) override;
  CheckerType* checker_ptr() const { return checker_ptr_.get(); }

 protected:
  /*! \brief Do the standard checks for convergence.
   *
   * This includes status of convergence, delta until convergence, and
   * iterates.
   *
   * \param current_iteration current iteration
   * \param previous_iteration previous iteration
   */
  void StatusDeltaAndIterate(CompareType& current_iteration,
                             CompareType& previous_iteration);

  using Final<CompareType>::convergence_status_;
  std::unique_ptr<CheckerType> checker_ptr_;
  std::shared_ptr<reporter::MpiI> reporter_ptr_ = nullptr;
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_MOMENTS_FINAL_OR_N_H_