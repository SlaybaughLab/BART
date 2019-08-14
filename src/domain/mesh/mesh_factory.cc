#include "domain/mesh/mesh_factory.h"

#include "domain/mesh/mesh_i.h"
#include "domain/mesh/mesh_cartesian.h"

namespace bart {

namespace domain {

namespace mesh {

template<int dim>
std::unique_ptr<MeshI<dim>> MeshFactory<dim>::MakeCartesianMesh(
    const std::vector<double> spatial_max,
    const std::vector<int> n_cells,
    const std::string material_mapping,
    CartesianMeshImpl type) {

  std::unique_ptr<MeshI<dim>> return_ptr = nullptr;

  if (type == CartesianMeshImpl::kDefault) {
    if (material_mapping != "") {
      return_ptr = std::move(
          std::make_unique<MeshCartesian<dim>>(
              spatial_max, n_cells, material_mapping));
    } else {
      return_ptr = std::move(
          std::make_unique<MeshCartesian<dim>>(spatial_max, n_cells));
    }

  }

  return return_ptr;
}

template class MeshFactory<1>;
template class MeshFactory<2>;
template class MeshFactory<3>;

} // namespace mesh

} // namespace domain

} // namespace bart