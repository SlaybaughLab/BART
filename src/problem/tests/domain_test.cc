#include "../domain.h"

#include <memory>

#include <gtest/gtest.h>

#include "../../test_helpers/gmock_wrapper.h"
#include "../../domain/tests/mesh_mock.h"
#include "../../domain/tests/finite_element_mock.h"

using ::testing::_;
using ::testing::NiceMock;

class DomainTest : public ::testing::Test {
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

void DomainTest::SetUp() {
  mesh_mock = std::make_unique<bart::domain::MeshMock<2>>();
  nice_mesh_mock = std::make_unique<NiceMock<bart::domain::MeshMock<2>>>();
  fe_mock = std::make_unique<bart::domain::FiniteElementMock<2>>();
}

TEST_F(DomainTest, Constructor) {
  EXPECT_CALL(*mesh_mock, FillTriangulation(_));
  EXPECT_CALL(*mesh_mock, FillMaterialID(_));
  EXPECT_CALL(*mesh_mock, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*mesh_mock, FillBoundaryID(_));

  MocksToPointers();
  
  bart::problem::Domain<2> test_domain(mesh_ptr, fe_ptr);
  EXPECT_EQ(mesh_ptr, nullptr);
  EXPECT_EQ(fe_ptr, nullptr);
}

TEST_F(DomainTest, ConstructorErrorMaterialMapping) {
  EXPECT_CALL(*nice_mesh_mock, has_material_mapping()).
      WillOnce(::testing::Return(false));

  MocksToPointers();
  
  EXPECT_ANY_THROW({
      bart::problem::Domain<2> test_domain(nice_mesh_ptr, fe_ptr);
    });
}
