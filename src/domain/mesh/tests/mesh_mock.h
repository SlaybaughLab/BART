#ifndef BART_SRC_DOMAIN_MESH_MOCK_H_
#define BART_SRC_DOMAIN_MESH_MOCK_H_

#include "domain/mesh/mesh_i.h"

#include <array>

#include <deal.II/grid/tria.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart {

namespace domain {

namespace mesh {

/* \brief Provides a dealii Mesh */

template <int dim>
class MeshMock : public MeshI<dim> {
 public:
  MOCK_METHOD1_T(FillTriangulation, void(dealii::Triangulation<dim>&));

  MOCK_METHOD1_T(FillMaterialID, void(dealii::Triangulation<dim>&));

  MOCK_METHOD1_T(FillBoundaryID, void(dealii::Triangulation<dim>&));

  MOCK_CONST_METHOD0_T(has_material_mapping, bool());

  MOCK_CONST_METHOD0_T(spatial_max, std::array<double, dim>());

  MOCK_CONST_METHOD0_T(n_cells,  std::array<int, dim>());

  MOCK_METHOD(std::string, description, (), (const, override));
};

} // namespace mesh

} // namespace domain

} // namespace bart 

#endif // BART_SRC_DOMAIN_MESH_MOCK_H_
