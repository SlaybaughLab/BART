#include "cartesian_mesh.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include "../problem/parameter_types.h"

namespace bart {

namespace domain {

template <int dim>
CartesianMesh<dim>::CartesianMesh(const std::vector<double> spatial_max,
                                  const std::vector<int> n_cells,
                                  const std::string material_mapping)
    : CartesianMesh(spatial_max, n_cells) {
  ParseMaterialMap(material_mapping);
}


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
    
  for (int j = 0; j < static_cast<int>(y_vector.size()); ++j ) {
      
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
void CartesianMesh<dim>::FillMaterialID(dealii::Triangulation<dim> &to_fill) {
  for (auto cell = to_fill.begin_active(); cell != to_fill.end(); ++cell) {
    if (cell->is_locally_owned()) {
      int material_id = GetMaterialID(cell->center());
      cell->set_material_id(material_id);
    }
  }
}

template <int dim>
void CartesianMesh<dim>::FillBoundaryID(dealii::Triangulation<dim> &to_fill) {
  using Boundary = bart::problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  double zero_tol = 1.0e-14;
  
  for (auto cell = to_fill.begin_active(); cell != to_fill.end(); ++cell) {
    if (cell->is_locally_owned()) {
      for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
        auto face = cell->face(face_id);
        dealii::Point<dim> face_center = face->center();

        switch (dim) {
          case 2: {
            if (std::fabs(face_center[1]) < zero_tol) {
              face->set_boundary_id(static_cast<int>(Boundary::kYMin));
              break;
            } else if (std::fabs(face_center[1] - spatial_max_[1]) < zero_tol) {
              face->set_boundary_id(static_cast<int>(Boundary::kYMax));
              break;
            }
            [[fallthrough]];
          }
          // Fall through to check x-direction
          case 1: {
            if (std::fabs(face_center[0]) < zero_tol) {
              face->set_boundary_id(static_cast<int>(Boundary::kXMin));
              break;
            } else if (std::fabs(face_center[0] - spatial_max_[0]) < zero_tol) {
              face->set_boundary_id(static_cast<int>(Boundary::kXMax));
              break;
            }
            break;
          }
          default: {
            AssertThrow(false,
                        dealii::ExcMessage("Unsupported number of dimensions in FillBoundaryID"));
          }
        }
      }
    }
  }
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

template <int dim>
void SetupTriangulation(dealii::Triangulation<dim> &to_setup,
                        CartesianMesh<dim> &mesh) {
  mesh.FillTriangulation(to_setup);
  mesh.FillBoundaryID(to_setup);
  if (mesh.has_material_mapping())
    mesh.FillMaterialID(to_setup);
}


template class CartesianMesh<1>;
template class CartesianMesh<2>;
template void SetupTriangulation<1>(dealii::Triangulation<1>&,
                                    CartesianMesh<1>&);
template void SetupTriangulation<2>(dealii::Triangulation<2>&,
                                    CartesianMesh<2>&);

} // namespace domain

} // namespace bart 
