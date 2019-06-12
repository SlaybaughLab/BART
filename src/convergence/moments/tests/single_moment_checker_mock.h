#include "convergence/moments/single_moment_checker_i.h"

#include <optional>

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace moments {

class SingleMomentCheckerMock : public SingleMomentCheckerI {
 public:
  MOCK_METHOD2(CheckIfConverged, bool(const system::moments::MomentVector&,
      const system::moments::MomentVector&));
  MOCK_CONST_METHOD0(is_converged, bool());
  MOCK_METHOD1(SetMaxDelta, void(const double to_set));
  MOCK_CONST_METHOD0(max_delta, double());
  MOCK_CONST_METHOD0(delta, std::optional<double>());
};

} // namespace moments

} // namespace convergence

} // namespace bart