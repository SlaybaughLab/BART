#ifndef BART_SRC_CALCULATOR_TWO_GRID_TESTS_SPECTRAL_SHAPE_MOCK_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_TESTS_SPECTRAL_SHAPE_MOCK_HPP_

#include "acceleration/two_grid/spectral_shape/spectral_shape_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::acceleration::two_grid::spectral_shape {

class SpectralShapeMock : public SpectralShapeI {
 public:
  using SpectralShapeI::DealiiMatrix;
  MOCK_METHOD(std::vector<double>, CalculateSpectralShape, (const DealiiMatrix&, const DealiiMatrix&), (override));
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_CALCULATOR_TWO_GRID_TESTS_SPECTRAL_SHAPE_MOCK_HPP_
