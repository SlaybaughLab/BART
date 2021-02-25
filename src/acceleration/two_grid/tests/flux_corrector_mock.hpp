#ifndef BART_SRC_ACCELERATION_TWO_GRID_TESTS_FLUX_CORRECTOR_MOCK_HPP_
#define BART_SRC_ACCELERATION_TWO_GRID_TESTS_FLUX_CORRECTOR_MOCK_HPP_

#include "acceleration/two_grid/flux_corrector_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::acceleration::two_grid {

class FluxCorrectorMock : public FluxCorrectorI {
 public:
  MOCK_METHOD(void, CorrectFlux, (Vector& flux_to_correct, const Vector& error, const int group), (const, override));
};

} // namespace bart::acceleration::two_grid

#endif //BART_SRC_ACCELERATION_TWO_GRID_TESTS_FLUX_CORRECTOR_MOCK_HPP_
