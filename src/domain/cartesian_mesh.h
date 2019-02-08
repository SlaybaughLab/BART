#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <vector>

#include <deal.II/grid/tria.h>

namespace bart {

namespace domain {

template <int dim>
class CartesianMesh {
 public:
  CartesianMesh() = default;
  ~CartesianMesh() = default;

  void FillTriangulation(dealii::Triangulation<dim> &to_fill,
                         std::vector<double> spatial_max,
                         std::vector<int> n_cells);
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
