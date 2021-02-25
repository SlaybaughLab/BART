#ifndef BART_SRC_DOMAIN_MESH_CARTESIAN_HPP_
#define BART_SRC_DOMAIN_MESH_CARTESIAN_HPP_

#include <array>
#include <string>
#include <map>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

#include "domain/mesh/mesh_i.hpp"

namespace bart::domain::mesh {

template <int dim>
class MeshCartesian : public MeshI<dim> {
 public:
  MeshCartesian(const std::vector<double> spatial_max, const std::vector<int> n_cells);
  MeshCartesian(const std::vector<double> spatial_max, const std::vector<int> n_cells,
                const std::string material_mapping);
  ~MeshCartesian() = default;

  auto FillTriangulation(dealii::Triangulation<dim> &to_fill) -> void override;
  /*! \brief Parses a material mapping string that indicates relative locations
   * of materials */
  auto ParseMaterialMap(std::string material_mapping) -> void ;
  auto FillMaterialID(dealii::Triangulation<dim> &to_fill) -> void override;
  auto FillBoundaryID(dealii::Triangulation<dim> &to_fill) -> void override;

  /*! \brief Get the material ID for a given location (array)*/
  auto GetMaterialID(std::array<double, dim> location) -> int;
    /*! \brief Get the material ID for a given location (dealii Point)*/
  auto GetMaterialID(dealii::Point<dim> location) -> int;
  auto has_material_mapping() const -> bool override { return !material_mapping_.empty(); };
  /*! \brief Get spatial maximum in each direction */
  auto spatial_max() const -> std::array<double, dim> override { return spatial_max_; };
  /*! \brief Get number of cells in each direction */
  auto n_cells() const -> std::array<int, dim> override { return n_cells_; };
 private:
  std::array<double, dim> spatial_max_;
  std::array<int, dim>    n_material_cells_;
  std::array<int, dim>    n_cells_;
  std::map<std::array<int, dim>, int> material_mapping_;
  static bool is_registered_;
};

} // namespace bart::domain::mesh

#endif // BART_SRC_DOMAIN_MESH_CARTESIAN_H_
