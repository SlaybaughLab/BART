#include "domain/mesh/mesh_factory.h"

#include "domain/mesh/mesh_cartesian.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Eq, ::testing::Pointwise;

template <typename DimensionWrapper>
class DomainMeshFactoryIntegrationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

};

TYPED_TEST_CASE(DomainMeshFactoryIntegrationTest, bart::testing::AllDimensions);

TYPED_TEST(DomainMeshFactoryIntegrationTest, MakeMesh) {
  constexpr int dim = this->dim;

  std::vector<double> spatial_max{btest::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 1, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};
  std::string material_mapping = "1 1";

  auto mesh_ptr =
      domain::mesh::MeshFactory<dim>::MakeCartesianMesh(spatial_max, n_cells, material_mapping);

  using ExpectedType = domain::mesh::MeshCartesian<dim>;

  ASSERT_NE(mesh_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(mesh_ptr.get()));
  EXPECT_THAT(mesh_ptr->spatial_max(), Pointwise(Eq(), spatial_max));
  EXPECT_THAT(mesh_ptr->n_cells(), Pointwise(Eq(), n_cells));
  EXPECT_TRUE(mesh_ptr->has_material_mapping());
}

TYPED_TEST(DomainMeshFactoryIntegrationTest, MakeMeshNoMapping) {
  constexpr int dim = this->dim;

  std::vector<double> spatial_max{btest::RandomVector(dim, 0, 100)};
  std::vector<double> n_cells_double{btest::RandomVector(dim, 1, 20)};
  std::vector<int> n_cells{n_cells_double.begin(), n_cells_double.end()};

  auto mesh_ptr =
      domain::mesh::MeshFactory<dim>::MakeCartesianMesh(spatial_max, n_cells);

  using ExpectedType = domain::mesh::MeshCartesian<dim>;

  ASSERT_NE(mesh_ptr, nullptr);
  ASSERT_NE(nullptr, dynamic_cast<ExpectedType*>(mesh_ptr.get()));
  EXPECT_THAT(mesh_ptr->spatial_max(), Pointwise(Eq(), spatial_max));
  EXPECT_THAT(mesh_ptr->n_cells(), Pointwise(Eq(), n_cells));
  EXPECT_FALSE(mesh_ptr->has_material_mapping());
}

} // namespace