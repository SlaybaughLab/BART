#include "cartesian_mesh.h"

#include <algorithm>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/tensor.h>

namespace bart {

namespace domain {

template <int dim>
CartesianMesh<dim>::CartesianMesh(const std::vector<double> spatial_max,
                                  const std::vector<int> n_cells) {

  // Check lengths of spatial max and n_cells
  AssertThrow(spatial_max.size() == dim,
              dealii::ExcMessage("CartesianMesh argument error, incorrect spatial vector size"));
  AssertThrow(n_cells.size() == dim,
              dealii::ExcMessage("CartesianMesh argument error, incorrect number of cells vector size"))
  
  std::copy(n_cells.begin(), n_cells.end(), n_cells_.begin());
  std::copy(spatial_max.begin(), spatial_max.end(), spatial_max_.begin());
}

template <int dim>
void CartesianMesh<dim>::FillTriangulation(dealii::Triangulation<dim> &to_fill) {
  
  dealii::Point<dim> diagonal;
  dealii::Point<dim> origin;

  for (int i = 0; i < dim; ++i)
    diagonal[i] = spatial_max_[i];

  std::vector<unsigned int> number_of_cells{n_cells_.begin(), n_cells_.end()};

  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}
template <int dim>
void CartesianMesh<dim>::SetMaterialIDs(dealii::Triangulation<dim> &to_set,
                                        std::string material_mapping) {};

template <int dim>
int CartesianMesh<dim>::GetMaterialID(std::array<double, dim> location) {return 0;};


template class CartesianMesh<1>;
template class CartesianMesh<2>;
template class CartesianMesh<3>;

} // namespace domain

} // namespace bart 
