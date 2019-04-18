#include "domain/definition.h"

#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>

#include "test_helpers/gmock_wrapper.h"
#include "domain/tests/mesh_mock.h"
#include "domain/tests/finite_element_mock.h"

namespace {

using ::testing::_;
using ::testing::NiceMock;

class DefinitionTest : public ::testing::Test {
 protected:

  std::unique_ptr<bart::domain::MeshMock<2>> mesh_ptr;
  std::unique_ptr<NiceMock<bart::domain::MeshMock<2>>> nice_mesh_ptr;
  std::shared_ptr<bart::domain::FiniteElementMock<2>> fe_ptr;

  void SetUp() override;
};

void DefinitionTest::SetUp() {
  mesh_ptr = std::make_unique<bart::domain::MeshMock<2>>();
  nice_mesh_ptr = std::make_unique<NiceMock<bart::domain::MeshMock<2>>>();
  fe_ptr = std::make_shared<bart::domain::FiniteElementMock<2>>();
}

TEST_F(DefinitionTest, Constructor) {
  bart::domain::Definition<2> test_domain(std::move(mesh_ptr), fe_ptr);
  // Verify ownership has been taken by constructor
  EXPECT_EQ(mesh_ptr, nullptr);
  EXPECT_EQ(fe_ptr.use_count(), 2);
}

TEST_F(DefinitionTest, SetUpMesh) {
  EXPECT_CALL(*mesh_ptr, FillTriangulation(_));
  EXPECT_CALL(*mesh_ptr, FillMaterialID(_));
  EXPECT_CALL(*mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*mesh_ptr, FillBoundaryID(_));

  bart::domain::Definition<2> test_domain(std::move(mesh_ptr), fe_ptr);

  EXPECT_NO_THROW(test_domain.SetUpMesh(););
}

TEST_F(DefinitionTest, SetUpMeshMaterialMappingError) {
  EXPECT_CALL(*nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(false));

  bart::domain::Definition<2> test_domain(std::move(nice_mesh_ptr), fe_ptr);
  EXPECT_ANY_THROW(test_domain.SetUpMesh(););
}

class DOFTest : public DefinitionTest {
 protected:
  DOFTest() : fe(1) {};
  dealii::Triangulation<2> triangulation;
  dealii::FE_Q<2> fe;
  int n_cells_;

  void SetUp() override;
  static void SetTriangulation(dealii::Triangulation<2> &to_fill) {
    dealii::GridGenerator::hyper_cube(to_fill, -1, 1);
    to_fill.refine_global(2);
  }
};

void DOFTest::SetUp() {
  DefinitionTest::SetUp();
}

TEST_F(DOFTest, SetUpDOFTest) {
  EXPECT_CALL(*nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*nice_mesh_ptr, FillTriangulation(_))
      .WillOnce(::testing::Invoke(SetTriangulation));
  EXPECT_CALL(*fe_ptr, finite_element())
      .WillOnce(::testing::Return(&fe));
  EXPECT_CALL(*fe_ptr, dofs_per_cell())
      .WillOnce(::testing::Return(4));

  bart::domain::Definition<2> test_domain(std::move(nice_mesh_ptr), fe_ptr);
  test_domain.SetUpMesh();
  test_domain.SetUpDOF();

  EXPECT_EQ(test_domain.total_degrees_of_freedom(), 25);
  EXPECT_EQ(test_domain.Cells().size(), 16);

  auto matrix = test_domain.GetCellMatrix();
  EXPECT_EQ(matrix.n_rows(), 4);
  EXPECT_EQ(matrix.n_cols(), 4);


}

} // namespace