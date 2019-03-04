#ifndef BART_SRC_CONVERGENCE_FLUX_TESTS_SINGLE_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_FLUX_TESTS_SINGLE_CHECKER_MOCK_H_

#include <optional>

#include "convergence/flux/single_checker_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace flux {

class SingleCheckerMock : public SingleCheckerI {
 public:
  MOCK_METHOD2(CheckIfConverged, bool(data::FluxVector&, data::FluxVector&));
  MOCK_CONST_METHOD0(is_converged, bool());
  MOCK_METHOD1(SetMaxDelta, void(double));
  MOCK_CONST_METHOD0(max_delta, double());
  MOCK_CONST_METHOD0(delta, std::optional<double>());
};

} // namespace flux

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_FLUX_TESTS_SINGLE_CHECKER_MOCK_H_
