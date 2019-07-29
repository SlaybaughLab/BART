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

template <typename DimensionWrapper>
class DefinitionTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  std::unique_ptr<bart::domain::MeshMock<dim>> mesh_ptr;
  std::unique_ptr<NiceMock<bart::domain::MeshMock<dim>>> nice_mesh_ptr;
  std::shared_ptr<bart::domain::FiniteElementMock<dim>> fe_ptr;

  void SetUp() override;
};

TYPED_TEST_CASE(DefinitionTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void DefinitionTest<DimensionWrapper>::SetUp() {
  mesh_ptr = std::make_unique<bart::domain::MeshMock<dim>>();
  nice_mesh_ptr = std::make_unique<NiceMock<bart::domain::MeshMock<dim>>>();
  fe_ptr = std::make_shared<bart::domain::FiniteElementMock<dim>>();
}

TYPED_TEST(DefinitionTest, Constructor) {
  bart::domain::Definition<this->dim> test_domain(std::move(this->mesh_ptr),
                                                  this->fe_ptr);
  // Verify ownership has been taken by constructor
  EXPECT_EQ(this->mesh_ptr, nullptr);
  EXPECT_EQ(this->fe_ptr.use_count(), 2);
}

TYPED_TEST(DefinitionTest, SetUpMesh) {
  EXPECT_CALL(*this->mesh_ptr, FillTriangulation(_));
  EXPECT_CALL(*this->mesh_ptr, FillMaterialID(_));
  EXPECT_CALL(*this->mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*this->mesh_ptr, FillBoundaryID(_));

  bart::domain::Definition<this->dim> test_domain(std::move(this->mesh_ptr),
                                          this->fe_ptr);

  EXPECT_NO_THROW(test_domain.SetUpMesh(););
}

TYPED_TEST(DefinitionTest, SetUpMeshMaterialMappingError) {
  EXPECT_CALL(*this->nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(false));

  bart::domain::Definition<this->dim> test_domain(std::move(this->nice_mesh_ptr),
                                            this->fe_ptr);
  EXPECT_ANY_THROW(test_domain.SetUpMesh(););
}

template <typename DimensionWrapper>
class DOFTest : public DefinitionTest<DimensionWrapper> {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  DOFTest() : fe(1) {};
  dealii::Triangulation<dim> triangulation;
  dealii::FE_Q<dim> fe;
  int n_cells_;

  void SetUp() override;
  static void SetTriangulation(dealii::Triangulation<dim> &to_fill) {
    dealii::GridGenerator::hyper_cube(to_fill, -1, 1);
    to_fill.refine_global(2);
  }
};

TYPED_TEST_CASE(DOFTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void DOFTest<DimensionWrapper>::SetUp() {
  DefinitionTest<DimensionWrapper>::SetUp();
}

TYPED_TEST(DOFTest, SetUpDOFTest) {
  EXPECT_CALL(*this->nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*this->nice_mesh_ptr, FillTriangulation(_))
      .WillOnce(::testing::Invoke(this->SetTriangulation));
  EXPECT_CALL(*this->fe_ptr, finite_element())
      .WillOnce(::testing::Return(&this->fe));
  EXPECT_CALL(*this->fe_ptr, dofs_per_cell())
      .Times(2)
      .WillRepeatedly(::testing::Return(4));

  bart::domain::Definition<this->dim> test_domain(std::move(this->nice_mesh_ptr),
                                                  this->fe_ptr);
  test_domain.SetUpMesh();
  test_domain.SetUpDOF();

  EXPECT_EQ(test_domain.total_degrees_of_freedom(), 25);
  EXPECT_EQ(test_domain.Cells().size(), 16);

  auto matrix = test_domain.GetCellMatrix();
  EXPECT_EQ(matrix.n_rows(), 4);
  EXPECT_EQ(matrix.n_cols(), 4);

  auto vector = test_domain.GetCellVector();
  EXPECT_EQ(vector.size(), 4);


}

} // namespace