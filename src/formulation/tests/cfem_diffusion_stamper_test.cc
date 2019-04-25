#include "formulation/cfem_diffusion_stamper.h"

#include <map>
#include <memory>
#include <type_traits>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "domain/tests/definition_mock.h"
#include "formulation/scalar/tests/cfem_diffusion_mock.h"
#include "problem/parameter_types.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"
#include "test_helpers/test_helper_functions.h"
#include "test_helpers/dealii_test_domain.h"

namespace {

using ::testing::AssertionResult;
using ::testing::AssertionFailure, ::testing::AssertionSuccess;
using ::testing::DoDefault;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::UnorderedElementsAreArray;
using ::testing::Unused;
using ::testing::_;

using namespace bart;
using bart::testing::CompareMPIMatrices;
using Cell = domain::DefinitionI<2>::Cell;
using InitToken = formulation::scalar::CFEM_DiffusionI<2>::InitializationToken;
using Matrix = dealii::FullMatrix<double>;

// TODO(Josh) Move this to a header where other stamper tests can use it
void OnesFill(Matrix& to_fill) {
  for (unsigned int i = 0; i < to_fill.n_rows(); ++i) {
    for (unsigned int j = 0; j < to_fill.n_cols(); ++j) {
      to_fill(i,j) += 1;
    }
  }
}

void FillMatrixWithOnes(Matrix& to_fill, Unused, Unused, Unused) {
  OnesFill(to_fill);
}

void FillMatrixWithOnesBoundary(Matrix& to_fill, Unused, Unused, Unused, Unused) {
  OnesFill(to_fill);
}

template <typename DimensionWrapper>
class CFEMDiffusionStamperTestExperimental : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value;
  using Cell = typename domain::DefinitionI<dim>::Cell;
  using InitToken = typename
      formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;
  std::unique_ptr<NiceMock<domain::DefinitionMock<dim>>> mock_definition_ptr;
  std::unique_ptr<NiceMock<formulation::scalar::CFEM_DiffusionMock<dim>>> mock_diffusion_ptr;
  void SetUp() override;
  InitToken init_token_;
};

template <typename DimensionWrapper>
void CFEMDiffusionStamperTestExperimental<DimensionWrapper>::SetUp() {
  mock_definition_ptr = std::make_unique<NiceMock<domain::DefinitionMock<dim>>>();
  mock_diffusion_ptr =
      std::make_unique<NiceMock<formulation::scalar::CFEM_DiffusionMock<dim>>>();

  ON_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillByDefault(Return(init_token_));

  Cell test_cell;
  std::vector<Cell> cells{test_cell};

  ON_CALL(*mock_definition_ptr, Cells())
      .WillByDefault(Return(cells));
}


TYPED_TEST_CASE(CFEMDiffusionStamperTestExperimental,
    bart::testing::AllDimensions);

TYPED_TEST(CFEMDiffusionStamperTestExperimental, Constructor) {
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;
  auto& dim = this->dim;

  EXPECT_CALL(*mock_definition_ptr, Cells())
      .WillOnce(DoDefault());
  EXPECT_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  EXPECT_EQ(mock_diffusion_ptr, nullptr);
  EXPECT_EQ(mock_definition_ptr, nullptr);

  EXPECT_TRUE(test_stamper.reflective_boundaries().empty());
}

TYPED_TEST(CFEMDiffusionStamperTestExperimental, SetReflective) {
  using Boundary = bart::problem::Boundary;
  auto& dim = this->dim;
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;

  std::unordered_set<problem::Boundary> reflective_boundaries = {
      problem::Boundary::kXMin,
  };

  formulation::CFEM_DiffusionStamper<dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.AddReflectiveBoundary(Boundary::kXMin);
  EXPECT_THAT(test_stamper.reflective_boundaries(),
              UnorderedElementsAreArray(reflective_boundaries));

  reflective_boundaries.erase(Boundary::kXMin);
  reflective_boundaries.insert(Boundary::kXMax);

  test_stamper.AddReflectiveBoundary(Boundary::kXMax)
      .RemoveReflectiveBoundary(Boundary::kXMin);

  EXPECT_THAT(test_stamper.reflective_boundaries(),
              UnorderedElementsAreArray(reflective_boundaries));

}

TYPED_TEST(CFEMDiffusionStamperTestExperimental, ConstructorWithReflective) {
  std::unordered_set<problem::Boundary> reflective_boundaries = {
      problem::Boundary::kYMax,
      problem::Boundary::kXMin,
  };

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(this->mock_diffusion_ptr),
      std::move(this->mock_definition_ptr),
      reflective_boundaries);

  EXPECT_THAT(test_stamper.reflective_boundaries(),
              UnorderedElementsAreArray(reflective_boundaries));

}

TYPED_TEST(CFEMDiffusionStamperTestExperimental, ConstructorWithReflectiveMap) {
  std::vector<problem::Boundary> reflective_boundaries = {
      problem::Boundary::kYMax,
      problem::Boundary::kXMin,
  };

  std::map<problem::Boundary, bool> reflective_boundary_map = {
      {problem::Boundary::kXMin, true},
      {problem::Boundary::kXMax, false},
      {problem::Boundary::kYMin, false},
      {problem::Boundary::kYMax, true}
  };

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper_2(
      std::move(this->mock_diffusion_ptr),
      std::move(this->mock_definition_ptr),
      reflective_boundary_map);

  EXPECT_THAT(test_stamper_2.reflective_boundaries(),
              UnorderedElementsAreArray(reflective_boundaries));
}


class CFEMDiffusionStamperTest : public ::testing::Test {
 protected:
  std::unique_ptr<NiceMock<domain::DefinitionMock<2>>> mock_definition_ptr;
  std::unique_ptr<NiceMock<formulation::scalar::CFEM_DiffusionMock<2>>> mock_diffusion_ptr;
  void SetUp() override;
  InitToken init_token_;

};


void CFEMDiffusionStamperTest::SetUp() {
  mock_definition_ptr = std::make_unique<NiceMock<domain::DefinitionMock<2>>>();
  mock_diffusion_ptr =
      std::make_unique<NiceMock<formulation::scalar::CFEM_DiffusionMock<2>>>();

  ON_CALL(*mock_diffusion_ptr, Precalculate(_))
      .WillByDefault(Return(init_token_));

  Cell test_cell;
  std::vector<Cell> cells{test_cell};

  ON_CALL(*mock_definition_ptr, Cells())
      .WillByDefault(Return(cells));
}


template <typename T>
class CFEMDiffusionStamperMPIExperiment
    : public ::testing::Test ,
      public T
{};


TYPED_TEST_CASE(CFEMDiffusionStamperMPIExperiment,
                bart::testing::DealiiTestDomains);

TYPED_TEST(CFEMDiffusionStamperMPIExperiment, Experiment) {
  EXPECT_TRUE(true);
}

// TODO(Josh) Put this in it's own header file?
class CFEMDiffusionStamperMPITests : public CFEMDiffusionStamperTest,
                                     public bart::testing::DealiiTestDomain<2> {
 protected:

  dealii::PETScWrappers::MPI::SparseMatrix& system_matrix_ = matrix_1;
  dealii::PETScWrappers::MPI::SparseMatrix& index_hits_ = matrix_2;
  dealii::PETScWrappers::MPI::SparseMatrix& boundary_hits_ = matrix_3;

  void SetUp() override;
  void SetUpBoundaries();
};

void CFEMDiffusionStamperMPITests::SetUp() {
  CFEMDiffusionStamperTest::SetUp();
  SetUpDealii();

  for (const auto& cell : cells_) {
    int mat_id = btest::RandomDouble(0, 10);
    cell->set_material_id(mat_id);
  }

  ON_CALL(*mock_definition_ptr, Cells())
      .WillByDefault(Return(cells_));
  dealii::FullMatrix<double> cell_matrix(fe_.dofs_per_cell,
                                         fe_.dofs_per_cell);
  ON_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillByDefault(Return(cell_matrix));

  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (auto cell : cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      for (auto index_j : local_dof_indices) {
        index_hits_.add(index_i, index_j, 1);
      }
    }
  }
  index_hits_.compress(dealii::VectorOperation::add);
}

TEST_F(CFEMDiffusionStamperMPITests, StampStreaming) {

  int group_number = 1;

  for (auto const& cell : cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
        FillCellStreamingTerm(_, _, cell, group_number))
        .WillOnce(Invoke(FillMatrixWithOnes));
  }
  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampStreamingTerm(system_matrix_, group_number);

  EXPECT_TRUE(CompareMPIMatrices(system_matrix_, index_hits_));
}

TEST_F(CFEMDiffusionStamperMPITests, StampCollision) {

  int group_number = 1;

  for (auto const& cell : cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
                FillCellCollisionTerm(_, _, cell, group_number))
        .WillOnce(Invoke(FillMatrixWithOnes));
  }
  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampCollisionTerm(system_matrix_, group_number);

  EXPECT_TRUE(CompareMPIMatrices(system_matrix_, index_hits_));
}

class CFEMDiffusionStamperBoundaryMPITests
    : public CFEMDiffusionStamperMPITests {
 protected:

  void SetUp() override;
};

void CFEMDiffusionStamperBoundaryMPITests::SetUp() {
  CFEMDiffusionStamperMPITests::SetUp();
  SetUpBoundaries();

  int faces_per_cell = dealii::GeometryInfo<2>::faces_per_cell;
  std::vector<dealii::types::global_dof_index> local_dof_indices(fe_.dofs_per_cell);

  for (auto &cell : cells_) {
    if (cell->at_boundary()) {
      cell->get_dof_indices(local_dof_indices);
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          for (auto index_i : local_dof_indices) {
            for (auto index_j : local_dof_indices) {
              boundary_hits_.add(index_i, index_j, 1);
            }
          }
        }
      }
    }
  }
  boundary_hits_.compress(dealii::VectorOperation::add);
}

void CFEMDiffusionStamperMPITests::SetUpBoundaries() {
  using Boundary = bart::problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<2>::faces_per_cell;
  double zero_tol = 1.0e-14;

  for (auto &cell : cells_) {
    for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
      auto face = cell->face(face_id);
      dealii::Point<2> face_center = face->center();

      if (std::fabs(face_center[1]) < zero_tol) {
        face->set_boundary_id(static_cast<int>(Boundary::kYMin));
      } else if (std::fabs(face_center[1] - 1) < zero_tol) {
        face->set_boundary_id(static_cast<int>(Boundary::kYMax));
      } else if (std::fabs(face_center[0]) < zero_tol) {
        face->set_boundary_id(static_cast<int>(Boundary::kXMin));
      } else if (std::fabs(face_center[0] - 1) < zero_tol) {
        face->set_boundary_id(static_cast<int>(Boundary::kXMax));
      }
    }
  }
}

TEST_F(CFEMDiffusionStamperBoundaryMPITests, StampVacuumBoundaryTerm) {
  using BoundaryType = formulation::scalar::CFEM_DiffusionI<2>::BoundaryType;
  int faces_per_cell = dealii::GeometryInfo<2>::faces_per_cell;

  for (auto const& cell : cells_) {
    if (cell->at_boundary()) {
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          EXPECT_CALL(*mock_diffusion_ptr,
                      FillBoundaryTerm(_, _, cell, face, BoundaryType::kVacuum))
              .WillOnce(Invoke(FillMatrixWithOnesBoundary));
        }
      }
    }
  }

  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampBoundaryTerm(system_matrix_);

  EXPECT_TRUE(CompareMPIMatrices(boundary_hits_, system_matrix_));
}

TEST_F(CFEMDiffusionStamperBoundaryMPITests, StampVacuumReflectiveBoundaryTerm) {
  using BoundaryType = formulation::scalar::CFEM_DiffusionI<2>::BoundaryType;
  using Boundary = problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<2>::faces_per_cell;

  for (auto const& cell : cells_) {
    if (cell->at_boundary()) {
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          if (cell->face(face)->boundary_id() == static_cast<int>(Boundary::kXMin)) {
            EXPECT_CALL(*mock_diffusion_ptr,
                        FillBoundaryTerm(_, _, cell, face, BoundaryType::kReflective))
                .WillOnce(Invoke(FillMatrixWithOnesBoundary));
          } else {
            EXPECT_CALL(*mock_diffusion_ptr,
                        FillBoundaryTerm(_, _, cell, face, BoundaryType::kVacuum))
                .WillOnce(Invoke(FillMatrixWithOnesBoundary));
          }
        }
      }
    }
  }

  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.AddReflectiveBoundary(Boundary::kXMin);

  test_stamper.StampBoundaryTerm(system_matrix_);
}



TEST_F(CFEMDiffusionStamperBoundaryMPITests, StampVacuumAndStreaming) {
  ON_CALL(*mock_diffusion_ptr, FillCellStreamingTerm(_, _, _, _))
      .WillByDefault(Invoke(FillMatrixWithOnes));
  ON_CALL(*mock_diffusion_ptr, FillBoundaryTerm(_, _, _, _, _))
      .WillByDefault(Invoke(FillMatrixWithOnesBoundary));

  formulation::CFEM_DiffusionStamper<2> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampStreamingTerm(system_matrix_, 0);
  test_stamper.StampBoundaryTerm(system_matrix_);

  index_hits_.add(1, boundary_hits_);

  EXPECT_TRUE(CompareMPIMatrices(index_hits_, system_matrix_));

}




} // namespace