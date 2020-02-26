#include "formulation/stamper.h"

#include "domain/tests/definition_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"

namespace {

using namespace bart;

using ::testing::ReturnRef, ::testing::Return, ::testing::DoDefault;

/* ===== BASIC TESTS ===========================================================
 * These tests verify basic functionality of formulation::Stamper. */
template <typename DimensionWrapper>
class FormulationStamperTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainDefinitionType = domain::DefinitionMock<dim>;

  std::shared_ptr<DomainDefinitionType> domain_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationStamperTest<DimensionWrapper>::SetUp() {
  domain_ptr_ = std::make_shared<DomainDefinitionType>();
}

TYPED_TEST_SUITE(FormulationStamperTest, bart::testing::AllDimensions);

// Constructor should set dependency correct, getter should return correct value
TYPED_TEST(FormulationStamperTest, Constructor) {
  using StamperType = formulation::Stamper<this->dim>;
  std::shared_ptr<StamperType> stamper_ptr;
  EXPECT_NO_THROW({
    stamper_ptr = std::make_shared<StamperType>(this->domain_ptr_);
  });
  ASSERT_NE(stamper_ptr->domain_ptr(), nullptr);
  EXPECT_EQ(stamper_ptr->domain_ptr(), this->domain_ptr_.get());
}
// Constructor should throw if a nullptr dependency is passed
TYPED_TEST(FormulationStamperTest, BadDependency) {
  using StamperType = formulation::Stamper<this->dim>;
  std::shared_ptr<StamperType> stamper_ptr;
  EXPECT_ANY_THROW({
    stamper_ptr = std::make_shared<StamperType>(nullptr);
  });
}

/* ===== DEAL.II DOMAIN TESTS ==================================================
 * These tests verify operation of the stamper on a real dealii domain in both
 * serial and MPI */
template <typename DimensionWrapper>
class FormulationStamperTestDealiiDomain
    : public ::testing::Test,
      public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using DomainDefinitionType = domain::DefinitionMock<dim>;
  using StamperType = formulation::Stamper<dim>;

  // Test object
  std::unique_ptr<StamperType> test_stamper_ptr_;

  // Depdendency
  std::shared_ptr<DomainDefinitionType> domain_ptr_;

  // Test parameters
  system::MPISparseMatrix& system_matrix = this->matrix_1;
  system::MPISparseMatrix& expected_matrix = this->matrix_2;
  system::MPISparseMatrix& boundary_expected_matrix = this->matrix_3;
  system::MPIVector& system_vector = this->vector_1;
  system::MPIVector& expected_vector = this->vector_2;
  system::MPIVector& boundary_expected_vector = this->vector_3;
  std::function<void(formulation::FullMatrix&,
                     const domain::CellPtr<dim>&)> matrix_stamp_function;
  std::function<void(formulation::Vector&,
                     const domain::CellPtr<dim>&)> vector_stamp_function;
  std::function<void(formulation::FullMatrix&,
                     const domain::FaceIndex,
                     const domain::CellPtr<dim>&)> matrix_boundary_stamp_function;
  std::function<void(formulation::Vector&,
                     const domain::FaceIndex,
                     const domain::CellPtr<dim>&)> vector_boundary_stamp_function;

  void SetUp() override;
};

void SetMatrixToOne(formulation::FullMatrix& to_stamp) {
  for (int i = 0; i < static_cast<int>(to_stamp.n_rows()); ++i) {
    for (int j = 0; j < static_cast<int>(to_stamp.n_cols()); ++j) {
      to_stamp(i,j) = 1;
    }
  }
}

void SetVectorToOne(formulation::Vector& to_set) {
  for (int i = 0; i < static_cast<int>(to_set.size()); ++i) {
    to_set(i) = 1;
  }
}

template <typename DimensionWrapper>
void FormulationStamperTestDealiiDomain<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
  domain_ptr_ = std::make_shared<DomainDefinitionType>();
  test_stamper_ptr_ = std::make_unique<StamperType>(domain_ptr_);

  /* Set up expected result of matrix stamping, we will return a matrix of all
   * ones for each cell. These should be stamped onto the degrees of freedom
   * for each cell. */
  int cell_dofs = this->fe_.dofs_per_cell;
  formulation::FullMatrix ones_matrix(cell_dofs, cell_dofs);
  SetMatrixToOne(ones_matrix);
  formulation::Vector ones_vector(cell_dofs);
  SetVectorToOne(ones_vector);

  for (const auto& cell : this->cells_) {
    std::vector<dealii::types::global_dof_index> local_dof_indices(cell_dofs);
    cell->get_dof_indices(local_dof_indices);
    expected_matrix.add(local_dof_indices, local_dof_indices, ones_matrix);
    expected_vector.add(local_dof_indices, ones_vector);
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          boundary_expected_vector.add(local_dof_indices, ones_vector);
          boundary_expected_matrix.add(local_dof_indices, local_dof_indices,
                                       ones_matrix);
        }
      }
    }
  }
  expected_matrix.compress(dealii::VectorOperation::add);
  expected_vector.compress(dealii::VectorOperation::add);
  boundary_expected_matrix.compress(dealii::VectorOperation::add);
  boundary_expected_vector.compress(dealii::VectorOperation::add);

  // Set up the stamp function that just sets the provided matrix to all 1
  matrix_stamp_function = [](formulation::FullMatrix& to_stamp,
                             const domain::CellPtr<dim>&) -> void {
    SetMatrixToOne(to_stamp);
  };
  vector_stamp_function = [](formulation::Vector& to_stamp,
                             const domain::CellPtr<dim>&) -> void {
    SetVectorToOne(to_stamp);
  };
  vector_boundary_stamp_function = [](formulation::Vector& to_stamp,
                                      const domain::FaceIndex,
                                      const domain::CellPtr<dim>&) -> void {
    SetVectorToOne(to_stamp);
  };
  matrix_boundary_stamp_function = [](formulation::FullMatrix& to_stamp,
                                      const domain::FaceIndex,
                                      const domain::CellPtr<dim>&) -> void {
    SetMatrixToOne(to_stamp);
  };

  // Set up expected calls for the definition object
  ON_CALL(*domain_ptr_, Cells()).WillByDefault(Return(this->cells_));
  ON_CALL(*domain_ptr_, GetCellMatrix())
      .WillByDefault(Return(dealii::FullMatrix<double>(cell_dofs, cell_dofs)));
  ON_CALL(*domain_ptr_, GetCellVector())
      .WillByDefault(Return(dealii::Vector<double>(cell_dofs)));
}

TYPED_TEST_SUITE(FormulationStamperTestDealiiDomain,
                 bart::testing::AllDimensions);

TYPED_TEST(FormulationStamperTestDealiiDomain, StampMatrixMPI) {
  EXPECT_CALL(*this->domain_ptr_, GetCellMatrix()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    this->test_stamper_ptr_->StampMatrix(this->system_matrix,
                                         this->matrix_stamp_function);
  });
  EXPECT_TRUE(test_helpers::CompareMPIMatrices(this->system_matrix,
                                               this->expected_matrix));
}

TYPED_TEST(FormulationStamperTestDealiiDomain, StampVectorMPI) {
  EXPECT_CALL(*this->domain_ptr_, GetCellVector()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    this->test_stamper_ptr_->StampVector(this->system_vector,
                                         this->vector_stamp_function);
                  });
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->system_vector,
                                              this->expected_vector));
}

TYPED_TEST(FormulationStamperTestDealiiDomain, StampMatrixBoundaryMPI) {
  EXPECT_CALL(*this->domain_ptr_, GetCellMatrix()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    this->test_stamper_ptr_->StampBoundaryMatrix(this->system_matrix,
                                                 this->matrix_boundary_stamp_function);
                  });
  EXPECT_TRUE(test_helpers::CompareMPIMatrices(this->system_matrix,
                                               this->boundary_expected_matrix));
}

TYPED_TEST(FormulationStamperTestDealiiDomain, StampVectorBoundaryMPI) {
  EXPECT_CALL(*this->domain_ptr_, GetCellVector()).WillOnce(DoDefault());
  EXPECT_CALL(*this->domain_ptr_, Cells()).WillOnce(DoDefault());
  EXPECT_NO_THROW({
    this->test_stamper_ptr_->StampBoundaryVector(this->system_vector,
                                                 this->vector_boundary_stamp_function);
                  });
  EXPECT_TRUE(test_helpers::CompareMPIVectors(this->system_vector,
                                              this->boundary_expected_vector));
}

} // namespace
