#ifndef BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_

#include "calculator/two_grid/material_spectral_shapes_i.hpp"
#include "test_helpers/gmock_wrapper.h"

namespace bart::calculator::two_grid {

class MaterialSpectralShapesMock : public MaterialSpectralShapesI {
 public:
  MOCK_METHOD(void, CalculateMaterialSpectralShapes, (std::shared_ptr<CrossSections>), (override));
  MOCK_METHOD((std::unordered_map<int, std::vector<double>>), material_spectral_shapes, (), (const, override));
};

} // namespace bart::calculator::two_grid

#endif //BART_SRC_CALCULATOR_TWO_GRID_TESTS_MATERIAL_SPECTRAL_SHAPES_MOCK_HPP_
