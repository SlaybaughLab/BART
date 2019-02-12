#include "cartesian_mesh.h"

#include <algorithm>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

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
void CartesianMesh<dim>::ParseMaterialMap(std::string material_mapping) {
  using dealii::Utilities::split_string_list;
  using dealii::Utilities::string_to_int;
  using StringVector = std::vector<std::string>;
      
  StringVector y_vector = split_string_list(material_mapping, "\n");
  
  std::reverse(y_vector.begin(), y_vector.end());
    
  for (int j = 0 ; j < static_cast<int>(y_vector.size()); ++j ) {
      
    StringVector x_vector = split_string_list(y_vector[j], " ");
      
    for (int i = 0; i < static_cast<int>(x_vector.size()); ++i ) {
        
      std::array<int, 2> location{i, j};
      material_mapping_[location] = string_to_int(x_vector[i]);
    }
    n_material_cells_[0] = x_vector.size();
  }
  n_material_cells_[1] = y_vector.size();
}

template <int dim>
int CartesianMesh<dim>::GetMaterialID(dealii::Point<dim> location) {
  std::array<double, dim> array_location;
  for (int i = 0; i < dim; ++i)
    array_location[i] = location[i];
  return GetMaterialID(array_location);
}
  
template <int dim>
int CartesianMesh<dim>::GetMaterialID(std::array<double, dim> location) {
  std::array<int, 2> relative_location{0, 0};

  for (int i = 0; i < dim; ++i) {
    double cell_size = spatial_max_[i]/n_material_cells_[i];
    double cell_location = location[i] / cell_size;
    int cell_index = std::floor(cell_location);

    if (static_cast<double>(cell_index) == cell_location && cell_index != 0) {
      relative_location[i] = cell_index - 1;
    } else {
      relative_location[i] = cell_index;         
    }
  }
  
  return material_mapping_[relative_location];
}


template class CartesianMesh<1>;
template class CartesianMesh<2>;

} // namespace domain

} // namespace bart 
