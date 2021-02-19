#ifndef BART_SRC_CONVERGENCE_TESTS_CONVERGENCE_CHECKER_MOCK_HPP_
#define BART_SRC_CONVERGENCE_TESTS_CONVERGENCE_CHECKER_MOCK_HPP_

#include "convergence/convergence_checker_i.hpp"

#include "test_helpers/gmock_wrapper.h"

namespace bart::convergence {

template <typename CompareT, typename DeltaT = double>
class ConvergenceCheckerMock : public ConvergenceCheckerI<CompareT, DeltaT> {
 public:
  MOCK_METHOD(bool, IsConverged, (const CompareT& current_iteration, const CompareT& previous_iteration), (override));
  MOCK_METHOD(bool, is_converged, (), (const, override));
  MOCK_METHOD(void, SetMaxDelta, (const DeltaT& to_set), (override));
  MOCK_METHOD(DeltaT, max_delta, (), (const override));
  MOCK_METHOD(std::optional<DeltaT>, delta, (), (const override));
  MOCK_METHOD(std::optional<int>, failed_index, (), (const, override));
};

} // namespace bart::convergence

#endif //BART_SRC_CONVERGENCE_TESTS_CONVERGENCE_CHECKER_MOCK_HPP_
