#ifndef BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_

#include <memory>

#include "material_spectral_shapes_i.hpp"
#include "spectral_shape_i.hpp"
#include "utility/has_dependencies.h"

namespace bart::acceleration::two_grid::spectral_shape {

class MaterialSpectralShapes : public MaterialSpectralShapesI, public utility::HasDependencies {
 public:
  using SpectralShapeCalculator = acceleration::two_grid::spectral_shape::SpectralShapeI;
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

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_HPP_
