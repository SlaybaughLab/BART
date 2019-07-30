#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_

#include <memory>

#include "convergence/final_i.h"
#include "iteration/outer/outer_iteration_i.h"
#include "iteration/updater/source_updater_i.h"


namespace bart {

namespace iteration {

namespace outer {

template <typename ConvergenceType>
class OuterIteration : public OuterIterationI {
 public:
  using ConvergenceChecker = convergence::FinalI<ConvergenceType>;
  using SourceUpdater = iteration::updater::SourceUpdaterI;

  OuterIteration(
      std::unique_ptr<ConvergenceChecker> convergence_checker_ptr,
      const std::shared_ptr<SourceUpdater> &source_updater_ptr);
  virtual ~OuterIteration() = default;
  virtual void IterateToConvergence(system::System &system);

  ConvergenceChecker* convergence_checker_ptr() const {
    return convergence_checker_ptr_.get();
  }

  SourceUpdater* source_updater_ptr() const {
    return source_updater_ptr_.get();
  }

 protected:
  virtual convergence::Status CheckConvergence(system::System &system) = 0;

  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_ = nullptr;
  std::shared_ptr<SourceUpdater> source_updater_ptr_ = nullptr;
};

} // namespace outer

} // namespace iteration

} // namespace bart

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_H_
