#ifndef BART_SRC_DOMAIN_MESH_I_H_
#define BART_SRC_DOMAIN_MESH_I_H_

#include <array>

#include <deal.II/grid/tria.h>

namespace bart {

namespace domain {
//! Classes providing a spatial mesh
namespace mesh {

/* \brief Provides a dealii Mesh */

template <int dim>
class MeshI {
 public:
  virtual ~MeshI() = default;

  /* \brief Generates the mesh in the triangulation object */
  virtual void FillTriangulation(dealii::Triangulation<dim> &to_fill) = 0;
  /* \brief Add material IDs to each cell in the triangulation */
  virtual void FillMaterialID(dealii::Triangulation<dim> &to_fill) = 0;
  /* \brief Add boundary IDs to each boundary face */
  virtual void FillBoundaryID(dealii::Triangulation<dim> &to_fill) = 0;

  /* \brief Gets if the triangulation has a material mapping */
  virtual bool has_material_mapping() const = 0;
  /* \brief Gets maximum spatial dimensions of the mesh */
  virtual std::array<double, dim> spatial_max() const = 0;
  /* \brief Gets number of cells in each dimension */
  virtual std::array<int, dim> n_cells() const = 0;

  virtual std::string description() const = 0;
};

} // namespace mesh

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_MESH_I_H_
