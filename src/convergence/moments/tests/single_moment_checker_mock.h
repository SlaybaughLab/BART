#include "convergence/moments/single_moment_checker_i.h"

#include <optional>

#include "test_helpers/gmock_wrapper.h"

namespace bart::convergence::moments {

class SingleMomentCheckerMock : public SingleMomentCheckerI {
 public:
  using MomentVector = system::moments::MomentVector;
  MOCK_METHOD(bool, IsConverged, (const MomentVector&, const MomentVector&), (override));
  MOCK_METHOD(bool, is_converged, (), (const, override));
  MOCK_METHOD(void, SetMaxDelta, (const double&), (override));
  MOCK_METHOD(double, max_delta, (), (const, override));
  MOCK_METHOD(std::optional<double>, delta, (), (const, override));
};

} // namespace bart::convergence::moments