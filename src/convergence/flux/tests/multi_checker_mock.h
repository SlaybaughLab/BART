#ifndef BART_SRC_CONVERGENCE_FLUX_TESTS_MULTI_CHECKER_MOCK_H_
#define BART_SRC_CONVERGENCE_FLUX_TESTS_MULTI_CHECKER_MOCK_H_

#include "convergence/flux/multi_checker_i.h"
#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

namespace flux {

class MultiCheckerMock : public MultiCheckerI {
 public:
  MOCK_METHOD2(CheckIfConverged, bool(data::MultiFluxPtrs &current_iteration,
                                      data::MultiFluxPtrs &previous_iteration));
  MOCK_CONST_METHOD0(is_converged, bool());
  MOCK_CONST_METHOD0(failed_index, std::optional<int>());
  MOCK_CONST_METHOD0(failed_delta, std::optional<double>());
};


} // namespace flux

} // namespace convergence

} // namespace bart


#endif // BART_SRC_CONVERGENCE_FLUX_TESTS_MULTI_CHECKER_MOCK_H_
