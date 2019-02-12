#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <array>
#include <string>
#include <map>
#include <vector>

#include <deal.II/base/point.h>
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
  void ParseMaterialMap(std::string material_mapping);
  void FillMaterialID(dealii::Triangulation<dim> &to_fill);
  void FillBoundaryID(dealii::Triangulation<dim> &to_fill);
  
  int GetMaterialID(std::array<double, dim> location);
  int GetMaterialID(dealii::Point<dim> location);
      
 private:
  std::array<double, dim> spatial_max_;
  std::array<int, 2>    n_material_cells_;
  std::array<int, dim>    n_cells_;
  std::map<std::array<int, 2>, int> material_mapping_;
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
