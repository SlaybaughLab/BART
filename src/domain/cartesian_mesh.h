#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <array>

#include <deal.II/grid/tria.h>

namespace bart {

namespace domain {

template <int dim>
class CartesianMesh {
 public:
  CartesianMesh() = default;
  ~CartesianMesh() = default;

  void FillTriangulation(dealii::Triangulation<dim> &to_fill,
                         std::array<double, dim> x_max,
                         std::array<int, dim> n_cells);
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
