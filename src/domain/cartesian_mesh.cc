#include "cartesian_mesh.h"

#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/tensor.h>

namespace bart {

namespace domain {

template <int dim>
void CartesianMesh<dim>::FillTriangulation(dealii::Triangulation<dim> &to_fill,
                                           std::vector<double> spatial_max,
                                           std::vector<int> n_cells) {
  dealii::Point<dim> diagonal;
  dealii::Point<dim> origin;

  for (int i = 0; i < dim; ++i)
    diagonal[i] = spatial_max[i];

  std::vector<unsigned int> number_of_cells{n_cells.begin(), n_cells.end()};

  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}


template class CartesianMesh<1>;
template class CartesianMesh<2>;
template class CartesianMesh<3>;

} // namespace domain

} // namespace bart 
