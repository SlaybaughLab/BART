#include "domain/definition.h"

#include <cmath>
#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include <deal.II/grid/tria.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/grid/grid_generator.h>

#include "test_helpers/gmock_wrapper.h"
#include "domain/mesh/tests/mesh_mock.h"
#include "domain/finite_element/tests/finite_element_mock.h"
#include "problem/parameter_types.h"

namespace {

using namespace bart;

using ::testing::_;
using ::testing::NiceMock;

template <typename DimensionWrapper>
class DomainDefinitionTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  std::unique_ptr<bart::domain::mesh::MeshMock<dim>> mesh_ptr;
  std::unique_ptr<NiceMock<bart::domain::mesh::MeshMock<dim>>> nice_mesh_ptr;
  std::shared_ptr<bart::domain::finite_element::FiniteElementMock<dim>> fe_ptr;

  void SetUp() override;
};

TYPED_TEST_CASE(DomainDefinitionTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void DomainDefinitionTest<DimensionWrapper>::SetUp() {
  mesh_ptr = std::make_unique<bart::domain::mesh::MeshMock<dim>>();
  nice_mesh_ptr = std::make_unique<NiceMock<bart::domain::mesh::MeshMock<dim>>>();
  fe_ptr = std::make_shared<bart::domain::finite_element::FiniteElementMock<dim>>();
}

TYPED_TEST(DomainDefinitionTest, Constructor) {

  using Discretization = problem::DiscretizationType;
  std::array<Discretization, 2> discretizations{
    Discretization::kContinuousFEM, Discretization::kDiscontinuousFEM};

  for (const auto discretization : discretizations) {
    bart::domain::Definition<this->dim> test_domain(std::move(this->mesh_ptr),
                                                    this->fe_ptr,
                                                    discretization);
    // Verify ownership has been taken by constructor
    EXPECT_EQ(this->mesh_ptr, nullptr);
    EXPECT_EQ(this->fe_ptr.use_count(), 2);
    EXPECT_EQ(test_domain.discretization_type(), discretization);
  }
}

TYPED_TEST(DomainDefinitionTest, SetUpMesh) {
  EXPECT_CALL(*this->mesh_ptr, FillTriangulation(_));
  EXPECT_CALL(*this->mesh_ptr, FillMaterialID(_));
  EXPECT_CALL(*this->mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(true));
  EXPECT_CALL(*this->mesh_ptr, FillBoundaryID(_));

  bart::domain::Definition<this->dim> test_domain(std::move(this->mesh_ptr),
                                          this->fe_ptr);

  EXPECT_NO_THROW(test_domain.SetUpMesh(););
}

TYPED_TEST(DomainDefinitionTest, SetUpMeshMaterialMappingError) {
  EXPECT_CALL(*this->nice_mesh_ptr, has_material_mapping()).
      WillOnce(::testing::Return(false));

  bart::domain::Definition<this->dim> test_domain(std::move(this->nice_mesh_ptr),
                                            this->fe_ptr);
  EXPECT_ANY_THROW(test_domain.SetUpMesh(););
}

template <typename DimensionWrapper>
class DomainDefinitionDOFTest : public DomainDefinitionTest<DimensionWrapper> {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  DomainDefinitionDOFTest() : fe(1) {};
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

TYPED_TEST_CASE(DomainDefinitionDOFTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void DomainDefinitionDOFTest<DimensionWrapper>::SetUp() {
  DomainDefinitionTest<DimensionWrapper>::SetUp();
}

TYPED_TEST(DomainDefinitionDOFTest, SetUpDOFTestMPI) {
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

TYPED_TEST(DomainDefinitionDOFTest, SystemMatrixMPI) {
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

  auto system_matrix_ptr = test_domain.MakeSystemMatrix();

  ASSERT_NE(system_matrix_ptr, nullptr);
  EXPECT_EQ(system_matrix_ptr->n(), test_domain.locally_owned_dofs().size());
  EXPECT_EQ(system_matrix_ptr->m(), test_domain.locally_owned_dofs().size());
}

TYPED_TEST(DomainDefinitionDOFTest, SystemVectorMPI) {
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

  auto system_vector_ptr = test_domain.MakeSystemVector();

  ASSERT_NE(system_vector_ptr, nullptr);
  EXPECT_EQ(system_vector_ptr->size(), test_domain.locally_owned_dofs().size());
}

} // namespace