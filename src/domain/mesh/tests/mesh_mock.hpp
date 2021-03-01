#ifndef BART_SRC_DOMAIN_MESH_MOCK_HPP_
#define BART_SRC_DOMAIN_MESH_MOCK_HPP_

#include "domain/mesh/mesh_i.hpp"

#include <array>

#include <deal.II/grid/tria.h>

#include "test_helpers/gmock_wrapper.h"

namespace bart::domain::mesh {

template <int dim>
class MeshMock : public MeshI<dim> {
 public:
  MOCK_METHOD(void, FillTriangulation, (dealii::Triangulation<dim>&), (override));

  MOCK_METHOD(void, FillMaterialID, (dealii::Triangulation<dim>&), (override));

  MOCK_METHOD(void, FillBoundaryID, (dealii::Triangulation<dim>&), (override));

  MOCK_METHOD(bool, has_material_mapping, (), (override, const));

  MOCK_METHOD((std::array<double, dim>), spatial_max, (), (override, const));

  MOCK_METHOD((std::array<int, dim>), n_cells, (), (override, const));
};

} // namespace bart::domain::mesh

#endif // BART_SRC_DOMAIN_MESH_MOCK_HPP_
