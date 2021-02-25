#include "acceleration/two_grid/spectral_shape/domain_spectral_shapes.hpp"

namespace bart::acceleration::two_grid::spectral_shape {

template<int dim>
auto DomainSpectralShapes<dim>::CalculateDomainSpectralShapes(
    const MaterialToGroupSpectralShapeMap& material_to_group_spectral_shape_map,
    const Domain& domain) const -> GroupToDomainSpectralShapeMap {
  const auto cells{ domain.Cells() };
  const int total_dofs { domain.total_degrees_of_freedom() };
  auto cell_vector = domain.GetCellVector();
  std::vector<unsigned int> cell_dof_indices(cell_vector.size());
  GroupToDomainSpectralShapeMap return_mapping;

  for (auto& cell : cells) {
    const int material_id = cell->material_id();
    auto spectral_shape_by_group = material_to_group_spectral_shape_map.at(material_id);
    cell->get_dof_indices(cell_dof_indices);
    for (auto group = 0; group < spectral_shape_by_group.size(); ++group) {
      cell_vector = 0;
      cell_vector.add(spectral_shape_by_group.at(group));
      if (return_mapping.find(group) == return_mapping.end()) {
        return_mapping[group] = dealii::Vector<double>(total_dofs);
      }
      return_mapping[group].add(cell_dof_indices, cell_vector);
    }
  }

  return return_mapping;
}

template class DomainSpectralShapes<1>;
template class DomainSpectralShapes<2>;
template class DomainSpectralShapes<3>;

} // namespace bart::acceleration::two_grid::spectral_shape
