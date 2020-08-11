#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_

#include <memory>

#include "convergence/final_i.h"
#include "instrumentation/port.h"
#include "iteration/group/group_solve_iteration_i.h"
#include "iteration/outer/outer_iteration_i.h"
#include "system/system.h"

namespace bart {

namespace iteration {

namespace outer {

namespace data_names {
struct GroupConvergenceStatus;
struct Status;
struct IterationError;
using ConvergenceStatusPort = instrumentation::Port<convergence::Status, GroupConvergenceStatus>;
using StatusPort = instrumentation::Port<std::string, Status>;
using IterationErrorPort = instrumentation::Port<std::pair<int, double>, IterationError>;
} // namespace data_names


template <typename ConvergenceType>
class OuterIteration : public OuterIterationI,
                       public data_names::ConvergenceStatusPort,
                       public data_names::StatusPort,
                       public data_names::IterationErrorPort {
 public:
  using GroupIterator = iteration::group::GroupSolveIterationI;
  using ConvergenceChecker = convergence::FinalI<ConvergenceType>;

  OuterIteration(
      std::unique_ptr<GroupIterator> group_iterator_ptr,
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr);
  virtual ~OuterIteration() = default;
  virtual void IterateToConvergence(system::System &system);

  GroupIterator* group_iterator_ptr() const {
    return group_iterator_ptr_.get();
  }

  ConvergenceChecker* convergence_checker_ptr() const {
    return convergence_checker_ptr_.get();
  }

  std::vector<double> iteration_error() const override {
    return iteration_error_;
  }

 protected:
  virtual void InnerIterationToConvergence(system::System &system);
  virtual convergence::Status CheckConvergence(system::System &system) = 0;
  virtual void UpdateSystem(system::System& system, const int group,
                            const int angle) = 0;

  std::unique_ptr<GroupIterator> group_iterator_ptr_ = nullptr;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_ = nullptr;
  std::vector<double> iteration_error_{};
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
