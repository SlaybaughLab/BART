#ifndef BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_I_HPP_
#define BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_I_HPP_

#include "data/cross_sections/cross_sections_i.hpp"

namespace bart::acceleration::two_grid::spectral_shape {

/*! \brief Interface for classes that calculates and stores the spectral shape functions for all materials. */
class MaterialSpectralShapesI {
 public:
  using CrossSections = data::cross_sections::CrossSectionsI;
  virtual ~MaterialSpectralShapesI() = default;
  virtual auto CalculateMaterialSpectralShapes(std::shared_ptr<CrossSections>) -> void = 0;
  virtual auto material_spectral_shapes() const -> std::unordered_map<int, std::vector<double>> = 0;
};

} // namespace bart::acceleration::two_grid::spectral_shape

#endif //BART_SRC_CALCULATOR_TWO_GRID_MATERIAL_SPECTRAL_SHAPES_I_HPP_
