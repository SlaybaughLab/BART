#include "cartesian_mesh.h"

#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>

namespace bart {

namespace domain {

void CartesianMesh::FillTriangulation(dealii::Triangulation<1> &to_fill,
                                      std::array<double, 1> x_max,
                                      std::array<int, 1> n_cells) {
  dealii::Point<1> origin;
  dealii::Point<1> diagonal{x_max[0]};

  std::vector<unsigned int> number_of_cells{n_cells.begin(), n_cells.end()};

  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}

void CartesianMesh::FillTriangulation(dealii::Triangulation<2> &to_fill,
                                      std::array<double, 2> max_dim,
                                      std::array<int, 2> n_cells) {
  dealii::Point<2> origin;
  dealii::Point<2> diagonal{max_dim[0], max_dim[1]};

  std::vector<unsigned int> number_of_cells{n_cells.begin(), n_cells.end()};

  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}

void CartesianMesh::FillTriangulation(dealii::Triangulation<3> &to_fill,
                                      std::array<double, 3> max_dim,
                                      std::array<int, 3> n_cells) {
  dealii::Point<3> origin;
  dealii::Point<3> diagonal{max_dim[0], max_dim[1], max_dim[2]};

  std::vector<unsigned int> number_of_cells{n_cells.begin(), n_cells.end()};

  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}


} // namespace domain

} // namespace bart 
