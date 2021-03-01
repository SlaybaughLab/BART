#ifndef BART_SRC_DOMAIN_MESH_I_HPP_
#define BART_SRC_DOMAIN_MESH_I_HPP_

#include <array>

#include <deal.II/grid/tria.h>

#include "utility/has_description.h"

//! Classes providing a spatial mesh
namespace bart::domain::mesh {

/* \brief Provides a dealii Mesh */

template <int dim>
class MeshI : public utility::HasDescription {
 public:
  virtual ~MeshI() = default;

  /* \brief Generates the mesh in the triangulation object */
  virtual auto FillTriangulation(dealii::Triangulation<dim> &to_fill) -> void = 0;
  /* \brief Add material IDs to each cell in the triangulation */
  virtual auto FillMaterialID(dealii::Triangulation<dim> &to_fill) -> void = 0;
  /* \brief Add boundary IDs to each boundary face */
  virtual auto FillBoundaryID(dealii::Triangulation<dim> &to_fill) -> void = 0;

  /* \brief Gets if the triangulation has a material mapping */
  virtual auto has_material_mapping() const -> bool = 0;
  /* \brief Gets maximum spatial dimensions of the mesh */
  virtual auto spatial_max() const -> std::array<double, dim> = 0;
  /* \brief Gets number of cells in each dimension */
  virtual auto n_cells() const -> std::array<int, dim> = 0;
};

} // namespace bart::domain::mesh

#endif // BART_SRC_DOMAIN_MESH_I_HPP_
