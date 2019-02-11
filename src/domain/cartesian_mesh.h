#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <array>
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
      
 private:
  std::array<double, dim> spatial_max_;
  std::array<int, dim>    n_cells_;
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
