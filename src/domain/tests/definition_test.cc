#include "domain/definition.h"

#include <cmath>
#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>

#include "data/matrix_parameters.h"
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
    if (dim == 1) {
      to_fill.refine_global(4);
    } else {
      to_fill.refine_global(2);
    }
  }
};

TYPED_TEST_CASE(DOFTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void DOFTest<DimensionWrapper>::SetUp() {
  DefinitionTest<DimensionWrapper>::SetUp();
}

TYPED_TEST(DOFTest, SetUpDOFTestMPI) {
  EXPECT_CALL(*this->nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*this->nice_mesh_ptr, FillTriangulation(_))
      .WillOnce(::testing::Invoke(this->SetTriangulation));
  EXPECT_CALL(*this->fe_ptr, finite_element())
      .WillOnce(::testing::Return(&this->fe));

  bart::domain::Definition<this->dim> test_domain(std::move(this->nice_mesh_ptr),
                                                  this->fe_ptr);
  test_domain.SetUpMesh();
  test_domain.SetUpDOF();

  EXPECT_EQ(test_domain.total_degrees_of_freedom(),
            test_domain.dof_handler().n_dofs());

  int total_cells = 0;
  for (auto cell = test_domain.dof_handler().begin_active();
       cell != test_domain.dof_handler().end(); ++cell) {
    if (cell->is_locally_owned())
      ++total_cells;
  }

  EXPECT_EQ(test_domain.Cells().size(), total_cells);

  EXPECT_CALL(*this->fe_ptr, dofs_per_cell())
      .Times(2)
      .WillRepeatedly(::testing::Return(4));

  auto matrix = test_domain.GetCellMatrix();
  EXPECT_EQ(matrix.n_rows(), 4);
  EXPECT_EQ(matrix.n_cols(), 4);

  auto vector = test_domain.GetCellVector();
  EXPECT_EQ(vector.size(), 4);
}

TYPED_TEST(DOFTest, MatrixParametersMPI) {
  EXPECT_CALL(*this->nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*this->nice_mesh_ptr, FillTriangulation(_))
      .WillOnce(::testing::Invoke(this->SetTriangulation));
  EXPECT_CALL(*this->fe_ptr, finite_element())
      .WillOnce(::testing::Return(&this->fe));

  bart::domain::Definition<this->dim> test_domain(std::move(this->nice_mesh_ptr),
                                                  this->fe_ptr);
  test_domain.SetUpMesh();
  test_domain.SetUpDOF();

  bart::data::MatrixParameters test_parameters;
  test_domain.FillMatrixParameters(test_parameters,
      bart::problem::DiscretizationType::kContinuousFEM);

  EXPECT_EQ(test_parameters.rows, test_domain.locally_owned_dofs());
  EXPECT_EQ(test_parameters.columns, test_domain.locally_owned_dofs());
}

} // namespace