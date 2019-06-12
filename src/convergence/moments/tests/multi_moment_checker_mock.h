#ifndef BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_

#include <optional>

#include "convergence/moments/multi_moment_checker_i.h"
#include "system/moments/moment_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace moments {

class MultiMomentCheckerMock : public MultiMomentCheckerI {
 public:
  MOCK_METHOD2(CheckIfConverged, bool(const data::MomentsMap&,
      const data::MomentsMap&));
  MOCK_CONST_METHOD0(is_converged, bool());
  MOCK_CONST_METHOD0(failed_index, std::optional<int>());
  MOCK_CONST_METHOD0(delta, std::optional<double>());
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_