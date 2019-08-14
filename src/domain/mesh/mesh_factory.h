#ifndef BART_SRC_DOMAIN_MESH_MESH_FACTORY_H_
#define BART_SRC_DOMAIN_MESH_MESH_FACTORY_H_

#include <memory>
#include <string>
#include <vector>

namespace bart {

namespace domain {

namespace mesh {
template <int dim> class MeshI;

enum class CartesianMeshImpl {
  kDefault = 0,
};

template <int dim>
class MeshFactory {
 public:
  static std::unique_ptr<MeshI<dim>> MakeCartesianMesh(
      const std::vector<double> spatial_max,
      const std::vector<int> n_cells,
      const std::string material_mapping = "",
      CartesianMeshImpl type = CartesianMeshImpl::kDefault);
};

} // namespace mesh

} // namespace domain

} //namespace bart

#endif //BART_SRC_DOMAIN_MESH_MESH_FACTORY_H_
