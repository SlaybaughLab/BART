#ifndef BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_

#include <optional>

#include "convergence/moments/multi_moment_checker_i.h"
#include "system/moments/spherical_harmonic_types.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace moments {

class MultiMomentCheckerMock : public MultiMomentCheckerI {
 public:
  MOCK_METHOD2(IsConverged, bool(const system::moments::MomentsMap&,
      const system::moments::MomentsMap&));
  MOCK_CONST_METHOD0(is_converged, bool());
  MOCK_CONST_METHOD0(failed_index, std::optional<int>());
  MOCK_CONST_METHOD0(delta, std::optional<double>());
};

} // namespace moments

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_TESTS_MULTI_MOMENT_CHECKER_MOCK_H_