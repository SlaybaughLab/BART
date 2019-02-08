#ifndef BART_SRC_DOMAIN_CARTESIAN_MESH_H_
#define BART_SRC_DOMAIN_CARTESIAN_MESH_H_

#include <array>

#include <deal.II/grid/tria.h>

namespace bart {

namespace domain {

class CartesianMesh {
 public:
  CartesianMesh() = default;
  ~CartesianMesh() = default;

  void FillTriangulation(dealii::Triangulation<1> &to_fill,
                         std::array<double, 1> x_max,
                         std::array<int, 1> n_cells);

  void FillTriangulation(dealii::Triangulation<2> &to_fill,
                         std::array<double, 2> max_dim,
                         std::array<int, 2> n_cells);

  void FillTriangulation(dealii::Triangulation<3> &to_fill,
                         std::array<double, 3> max_dim,
                         std::array<int, 3> n_cells);
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_CARTESIAN_MESH_I_H_
