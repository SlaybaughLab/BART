#include "material_spectral_shapes.hpp"

#include <deal.II/lac/full_matrix.h>

namespace bart::acceleration::two_grid::spectral_shape {

MaterialSpectralShapes::MaterialSpectralShapes(std::unique_ptr<SpectralShapeCalculator> spectral_shape_calculator_ptr)
    : spectral_shape_calculator_ptr_(std::move(spectral_shape_calculator_ptr)) {
  AssertPointerNotNull(spectral_shape_calculator_ptr_.get(), "spectral shape calculator",
                       "MaterialSpectralShapesConstructor");
}
void MaterialSpectralShapes::CalculateMaterialSpectralShapes(std::shared_ptr<CrossSections> cross_sections_ptr) {
  using DealiiMatrix = dealii::FullMatrix<double>;
  const auto sigma_s_matrix_map{ cross_sections_ptr->sigma_s() };
  const auto sigma_t_vector_map{ cross_sections_ptr->sigma_t() };
  const int n_materials = sigma_t_vector_map.size();
  const int n_groups = sigma_t_vector_map.begin()->second.size();

  for (const auto& [material_id, sigma_t_vector] : sigma_t_vector_map) {
    DealiiMatrix sigma_t_matrix(n_groups, n_groups);
    for (int group = 0; group < n_groups; ++group) {
      sigma_t_matrix(group, group) = sigma_t_vector.at(group);
    }
    material_spectral_shapes_[material_id] =
        this->spectral_shape_calculator_ptr_->CalculateSpectralShape(sigma_t_matrix, sigma_s_matrix_map.at(material_id));
  }
}

} // namespace bart::acceleration::two_grid::spectral_shape
