#ifndef BART_SRC_DOMAIN_MESH_CARTESIAN_H_
#define BART_SRC_DOMAIN_MESH_CARTESIAN_H_

#include <array>
#include <string>
#include <map>
#include <vector>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

#include "mesh_i.h"

namespace bart {

namespace domain {

template <int dim>
class MeshCartesian : public MeshI<dim> {
 public:
  MeshCartesian(const std::vector<double> spatial_max,
                const std::vector<int> n_cells);
  MeshCartesian(const std::vector<double> spatial_max,
                const std::vector<int> n_cells,
                const std::string material_mapping);
  ~MeshCartesian() = default;

  void FillTriangulation(dealii::Triangulation<dim> &to_fill) override;
  /*! \brief Parses a material mapping string that indicates relative locations
   * of materials */
  void ParseMaterialMap(std::string material_mapping);
  void FillMaterialID(dealii::Triangulation<dim> &to_fill) override;
  void FillBoundaryID(dealii::Triangulation<dim> &to_fill) override;

  /*! \brief Get the material ID for a given location (array)*/
  int GetMaterialID(std::array<double, dim> location);
    /*! \brief Get the material ID for a given location (dealii Point)*/
  int GetMaterialID(dealii::Point<dim> location);
  bool has_material_mapping() const override {
    return !material_mapping_.empty(); };
  /*! \brief Get spatial maximum in each direction */
  std::array<double, dim> spatial_max() const override { return spatial_max_; };
  /*! \brief Get number of cells in each direction */
  std::array<int, dim> n_cells() const override { return n_cells_; };
      
 private:
  std::array<double, dim> spatial_max_;
  std::array<int, 2>    n_material_cells_;
  std::array<int, dim>    n_cells_;
  std::map<std::array<int, 2>, int> material_mapping_;
};

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_MESH_CARTESIAN_H_
