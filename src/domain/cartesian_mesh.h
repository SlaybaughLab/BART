#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <array>
#include <string>
#include <map>
#include <vector>

#include <deal.II/grid/tria.h>

namespace bart {

namespace domain {

template <int dim>
class CartesianMesh {
 public:
  CartesianMesh(const std::vector<double> spatial_max,
                const std::vector<int> n_cells);
  ~CartesianMesh() = default;

  void FillTriangulation(dealii::Triangulation<dim> &to_fill);
  void SetMaterialIDs(dealii::Triangulation<dim> &to_set,
                      std::string material_mapping);
  int GetMaterialID(std::array<double, dim> location);
      
 private:
  std::array<double, dim> spatial_max_;
  std::array<int, 3>    n_material_cells_;
  std::array<int, dim>    n_cells_;
  std::map<std::array<int, 3>, int> material_mapping_;
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
