#include "../definition.h"

#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>

#include "../../test_helpers/gmock_wrapper.h"
#include "../../domain/tests/mesh_mock.h"
#include "../../domain/tests/finite_element_mock.h"

using ::testing::_;
using ::testing::NiceMock;

class DefinitionTest : public ::testing::Test {
 protected:

  std::unique_ptr<bart::domain::MeshI<2>> mesh_ptr;
  std::unique_ptr<bart::domain::MeshI<2>> nice_mesh_ptr;
  std::unique_ptr<bart::domain::FiniteElementI<2>> fe_ptr;

  std::unique_ptr<bart::domain::MeshMock<2>> mesh_mock;
  std::unique_ptr<NiceMock<bart::domain::MeshMock<2>>> nice_mesh_mock;
  std::unique_ptr<bart::domain::FiniteElementMock<2>> fe_mock;
  
  void SetUp() override;
  void MocksToPointers() {
    mesh_ptr = std::move(mesh_mock);
    nice_mesh_ptr = std::move(nice_mesh_mock);
    fe_ptr = std::move(fe_mock); };
};

void DefinitionTest::SetUp() {
  mesh_mock = std::make_unique<bart::domain::MeshMock<2>>();
  nice_mesh_mock = std::make_unique<NiceMock<bart::domain::MeshMock<2>>>();
  fe_mock = std::make_unique<bart::domain::FiniteElementMock<2>>();
}

TEST_F(DefinitionTest, Constructor) {
  MocksToPointers();
  
  bart::domain::Definition<2> test_domain(mesh_ptr, fe_ptr);
  // Verify ownership has been taken by constructor
  EXPECT_EQ(mesh_ptr, nullptr);
  EXPECT_EQ(fe_ptr, nullptr);
}

TEST_F(DefinitionTest, SetUpMesh) {
  EXPECT_CALL(*mesh_mock, FillTriangulation(_));
  EXPECT_CALL(*mesh_mock, FillMaterialID(_));
  EXPECT_CALL(*mesh_mock, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*mesh_mock, FillBoundaryID(_));

  MocksToPointers();

  bart::domain::Definition<2> test_domain(mesh_ptr, fe_ptr);

  EXPECT_NO_THROW(test_domain.SetUpMesh(););
}
  
TEST_F(DefinitionTest, SetUpMeshMaterialMappingError) {
  EXPECT_CALL(*nice_mesh_mock, has_material_mapping()).
      WillOnce(::testing::Return(false));

  MocksToPointers();

  bart::domain::Definition<2> test_domain(nice_mesh_ptr, fe_ptr);
  EXPECT_ANY_THROW(test_domain.SetUpMesh(););
}

class DOFTest : public DefinitionTest {
 protected:
  DOFTest() : fe(1) {};
  dealii::Triangulation<2> triangulation;
  dealii::FE_Q<2>          fe;
  void SetUp() override;
  static void SetTriangulation(dealii::Triangulation<2> &to_fill) {
    dealii::GridGenerator::hyper_cube(to_fill, -1, 1);
    to_fill.refine_global(2);}
};

void DOFTest::SetUp() {
  DefinitionTest::SetUp();
}

TEST_F(DOFTest, SetUpDOFTest) {
  EXPECT_CALL(*nice_mesh_mock, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*nice_mesh_mock, FillTriangulation(_))
      .WillOnce(::testing::Invoke(SetTriangulation));
  EXPECT_CALL(*fe_mock, finite_element())
      .WillOnce(::testing::Return(&fe));

  MocksToPointers();
  bart::domain::Definition<2> test_domain(nice_mesh_ptr, fe_ptr);
  test_domain.SetUpMesh();
  test_domain.SetUpDOF();

  EXPECT_EQ(test_domain.total_degrees_of_freedom(), 25);
}
