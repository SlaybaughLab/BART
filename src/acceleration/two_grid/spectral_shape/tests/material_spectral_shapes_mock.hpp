#ifndef BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_

#include "acceleration/two_grid/spectral_shape/material_spectral_shapes_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::acceleration::two_grid::spectral_shape {

class MaterialSpectralShapesMock : public MaterialSpectralShapesI {
 public:
  MOCK_METHOD(void, CalculateMaterialSpectralShapes, (std::shared_ptr<CrossSections>), (override));
  MOCK_METHOD((std::unordered_map<int, std::vector<double>>), material_spectral_shapes, (), (const, override));
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_
