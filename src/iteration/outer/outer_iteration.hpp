#ifndef BART_SRC_ITERATION_OUTER_OUTER_ITERATION_HPP_
#define BART_SRC_ITERATION_OUTER_OUTER_ITERATION_HPP_

#include <memory>

#include "convergence/final_i.h"
#include "instrumentation/port.hpp"
#include "iteration/group/group_solve_iteration_i.h"
#include "iteration/outer/outer_iteration_i.hpp"
#include "system/system.hpp"
#include "system/moments/spherical_harmonic_i.h"
#include "utility/uncopyable.h"

namespace bart::iteration::outer {

namespace data_names {
using ConvergenceStatusPort = instrumentation::Port<convergence::Status, struct GroupConvergenceStatus>;
using StatusPort = instrumentation::Port<std::string, struct Status>;
using IterationErrorPort = instrumentation::Port<std::pair<int, double>, struct IterationError>;
using SolutionMomentsPort = instrumentation::Port<system::moments::SphericalHarmonicI, struct SolutionMoments>;
} // namespace data_names


template <typename ConvergenceType>
class OuterIteration : public OuterIterationI,
                       public utility::Uncopyable,
                       public data_names::ConvergenceStatusPort,
                       public data_names::StatusPort,
                       public data_names::IterationErrorPort,
                       public data_names::SolutionMomentsPort {
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

 protected:
  virtual void InnerIterationToConvergence(system::System &system);
  virtual convergence::Status CheckConvergence(system::System &system) = 0;
  virtual void UpdateSystem(system::System& system, const int group,
                            const int angle) = 0;

  std::unique_ptr<GroupIterator> group_iterator_ptr_ = nullptr;
  std::unique_ptr<ConvergenceChecker> convergence_checker_ptr_ = nullptr;
};

} // namespace bart::iteration::outer

#endif //BART_SRC_ITERATION_OUTER_OUTER_ITERATION_HPP_
