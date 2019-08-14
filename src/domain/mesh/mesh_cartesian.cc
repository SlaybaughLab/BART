#include "domain/mesh/mesh_cartesian.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include "problem/parameter_types.h"

namespace bart {

namespace domain {

namespace mesh {

template <int dim>
MeshCartesian<dim>::MeshCartesian(const std::vector<double> spatial_max,
                                  const std::vector<int> n_cells,
                                  const std::string material_mapping)
    : MeshCartesian(spatial_max, n_cells) {
  ParseMaterialMap(material_mapping);
}


template <int dim>
MeshCartesian<dim>::MeshCartesian(const std::vector<double> spatial_max,
                                  const std::vector<int> n_cells) {

  // Check lengths of spatial max and n_cells
  AssertThrow(spatial_max.size() == dim,
              dealii::ExcMessage("MeshCartesian argument error, incorrect spatial vector size"));
  AssertThrow(n_cells.size() == dim,
              dealii::ExcMessage("MeshCartesian argument error, incorrect number of cells vector size"));
  
  std::copy(n_cells.begin(), n_cells.end(), n_cells_.begin());
  std::copy(spatial_max.begin(), spatial_max.end(), spatial_max_.begin());
}

template <int dim>
void MeshCartesian<dim>::FillTriangulation(dealii::Triangulation<dim> &to_fill) {
  
  dealii::Point<dim> diagonal;
  dealii::Point<dim> origin;

  for (int i = 0; i < dim; ++i)
    diagonal[i] = spatial_max_[i];

  std::vector<unsigned int> number_of_cells{n_cells_.begin(), n_cells_.end()};
  dealii::GridGenerator::subdivided_hyper_rectangle(to_fill, number_of_cells,
                                                    origin, diagonal);
}

template <int dim>
void MeshCartesian<dim>::ParseMaterialMap(std::string material_mapping) {
  using dealii::Utilities::split_string_list;
  using dealii::Utilities::string_to_int;
  using StringVector = std::vector<std::string>;

  StringVector z_blocks = split_string_list(material_mapping, "\n\n");

  std::reverse(z_blocks.begin(), z_blocks.end());

  for (unsigned k = 0; k < z_blocks.size(); ++k) {
    StringVector y_line = split_string_list(z_blocks.at(k), "\n");
    std::reverse(y_line.begin(), y_line.end());
    for (unsigned j = 0; j < y_line.size(); ++j) {
      StringVector x_positions = split_string_list(y_line.at(j), " ");
      for (unsigned i = 0; i < x_positions.size(); ++i) {
        std::array<unsigned, 3> index{i, j, k};
        std::array<int, dim> location;
        for (int dir = 0; dir < dim; ++dir)
          location.at(dir) = index.at(dir);
        material_mapping_[location] = string_to_int(x_positions.at(i));
      }
      n_material_cells_.at(0) = x_positions.size();
    }
    if (dim > 1)
      n_material_cells_.at(1) = y_line.size();
  }

  if (dim > 2)
    n_material_cells_.at(2) = z_blocks.size();
}

template <int dim>  
void MeshCartesian<dim>::FillMaterialID(dealii::Triangulation<dim> &to_fill) {
  for (auto cell = to_fill.begin_active(); cell != to_fill.end(); ++cell) {
    if (cell->is_locally_owned()) {
      int material_id = GetMaterialID(cell->center());
      cell->set_material_id(material_id);
    }
  }
}

template <int dim>
void MeshCartesian<dim>::FillBoundaryID(dealii::Triangulation<dim> &to_fill) {
  using Boundary = bart::problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<dim>::faces_per_cell;
  double zero_tol = 1.0e-14;
  
  for (auto cell = to_fill.begin_active(); cell != to_fill.end(); ++cell) {
    if (cell->is_locally_owned() && cell->at_boundary()) {
      for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
        auto face = cell->face(face_id);
        if (face->at_boundary()) {
          dealii::Point<dim> face_center = face->center();
          switch (dim) {
            case 3: {
              if (std::fabs(face_center[2]) < zero_tol) {
                face->set_boundary_id(static_cast<int>(Boundary::kZMin));
                break;
              } else if (std::fabs(face_center[2] - spatial_max_.at(2)) < zero_tol) {
                face->set_boundary_id(static_cast<int>(Boundary::kZMax));
                break;
              }
              [[fallthrough]];
            }
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
}
  

template <int dim>
int MeshCartesian<dim>::GetMaterialID(dealii::Point<dim> location) {
  std::array<double, dim> array_location;
  for (int i = 0; i < dim; ++i)
    array_location[i] = location[i];
  return GetMaterialID(array_location);
}
  
template <int dim>
int MeshCartesian<dim>::GetMaterialID(std::array<double, dim> location) {
  std::array<int, dim> relative_location;

  for (int i = 0; i < dim; ++i) {
    double cell_size = spatial_max_.at(i)/n_material_cells_.at(i);
    double cell_location = location.at(i) / cell_size;
    int cell_index = std::floor(cell_location);

    if (static_cast<double>(cell_index) == cell_location && cell_index != 0) {
      relative_location.at(i) = cell_index - 1;
    } else {
      relative_location.at(i) = cell_index;
    }
  }
  
  return material_mapping_[relative_location];
}

template class MeshCartesian<1>;
template class MeshCartesian<2>;
template class MeshCartesian<3>;

} // namespace mesh

} // namespace domain

} // namespace bart 
