#ifndef BART_SRC_CONVERGENCE_TESTS_FLUX_MOCK_H_
#define BART_SRC_CONVERGENCE_TESTS_FLUX_MOCK_H_

#include "../flux_i.h"

#include "../../test_helpers/gmock_wrapper.h"

namespace bart {

namespace convergence {

class FluxMock : public FluxI {
 public:
  MOCK_METHOD2(isConverged, bool(data::Flux &, data::Flux &));
};

} // namespace convergence

} // namespace bart

#endif // BART_SRC_CONVERGENCE_TESTS_FLUX_MOCK_H_
