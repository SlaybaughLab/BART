#ifndef BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_

#include <memory>

#include "calculator/two_grid/material_spectral_shapes_i.hpp"
#include "calculator/two_grid/spectral_shape_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::calculator::two_grid {

class MaterialSpectralShapes : public MaterialSpectralShapesI, public utility::HasDependencies {
 public:
  using SpectralShapeCalculator = calculator::two_grid::SpectralShapeI;
  explicit MaterialSpectralShapes(std::unique_ptr<SpectralShapeCalculator>);

  auto CalculateMaterialSpectralShapes(std::shared_ptr<CrossSections> ptr) -> void override;
  auto material_spectral_shapes() const -> std::unordered_map<int, std::vector<double>> override {
    return material_spectral_shapes_;
  };

  auto spectral_shape_caluculator_ptr() { return spectral_shape_calculator_ptr_.get(); };
 private:
  std::unique_ptr<SpectralShapeCalculator> spectral_shape_calculator_ptr_{ nullptr };
  std::unordered_map<int, std::vector<double>> material_spectral_shapes_;
};

} // namespace bart::calculator::two_grid

#endif //BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_
