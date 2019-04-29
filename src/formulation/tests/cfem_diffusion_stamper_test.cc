#include "formulation/cfem_diffusion_stamper.h"

#include <map>
#include <memory>
#include <type_traits>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>

#include "data/moment_types.h"
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
using ::testing::Ref;
using ::testing::DoDefault;
using ::testing::Invoke;
using ::testing::NiceMock;
using ::testing::Return;
using ::testing::UnorderedElementsAreArray;
using ::testing::Unused;
using ::testing::_;

using namespace bart;
using bart::testing::CompareMPIMatrices;
using bart::testing::CompareMPIVectors;
using Matrix = dealii::FullMatrix<double>;
using Vector = dealii::Vector<double>;

// TODO(Josh) Move this to a header where other stamper tests can use it
void OnesFill(Matrix& to_fill) {
  for (unsigned int i = 0; i < to_fill.n_rows(); ++i) {
    for (unsigned int j = 0; j < to_fill.n_cols(); ++j) {
      to_fill(i,j) += 1;
    }
  }
}

void OnesFill(Vector& to_fill) {
  for (unsigned int i = 0; i < to_fill.size(); ++i) {
    to_fill(i) += 1;
  }
}

void FillMatrixWithOnes(Matrix& to_fill, Unused, Unused, Unused) {
  OnesFill(to_fill);
}

void FillMatrixWithOnesBoundary(Matrix& to_fill, Unused, Unused, Unused, Unused) {
  OnesFill(to_fill);
}

void FillVectorWithOnes3(Vector& to_fill, Unused, Unused) {
  OnesFill(to_fill);
}

void FillVectorWithOnes6(Vector& to_fill, Unused, Unused, Unused, Unused, Unused) {
  OnesFill(to_fill);
}

/* =============================================================================
 *
 * CFEM_DiffusionStamper Tests
 *
 * Tests the operation of the CFEM_DiffusionStamper. Test fixtures are divided
 * into the following:
 *
 * CFEMDiffusionStamperTest: contains pointers for the two dependencies of
 * CFEM_DiffusionStamper, but not much else. Templated to run in all three
 * dimensions.
 *
 * CFEMDiffusionStamperMPITests: derives from the previous, also sets up a full
 * dealii test domain to check that the stamper is stamping matrices and vectors
 * properly given MPI matrices or vectos and a list of cells.
 *
 * CFEMDiffusionStamperBoundaryMPITests: derives from the previous, checks the
 * stamping operation for boundary terms.
 *
 * =============================================================================
 */

/*
 * CFEMDiffusionStamperTest: Tests basic functionality
 *
 * Template type DimensionWrapper is intented to be a
 * std::integral_constant<int, N> where N is the dimension. This enables
 * automatic testing in all three dimensions.
 *
 */
template <typename DimensionWrapper>
class CFEMDiffusionStamperTest : public ::testing::Test {
 protected:
  static constexpr int dim = DimensionWrapper::value; // Problem dimension

  using Cell = typename domain::DefinitionI<dim>::Cell;
  using InitToken = typename
      formulation::scalar::CFEM_DiffusionI<dim>::InitializationToken;

  std::unique_ptr<NiceMock<domain::DefinitionMock<dim>>> mock_definition_ptr;
  std::unique_ptr<NiceMock<formulation::scalar::CFEM_DiffusionMock<dim>>> mock_diffusion_ptr;
  void SetUp() override;
  InitToken init_token_;
};

template <typename DimensionWrapper>
void CFEMDiffusionStamperTest<DimensionWrapper>::SetUp() {
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

// Initialize this as a typed test case, testing across 1D/2D/3D
TYPED_TEST_CASE(CFEMDiffusionStamperTest, bart::testing::AllDimensions);

TYPED_TEST(CFEMDiffusionStamperTest, Constructor) {
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

TYPED_TEST(CFEMDiffusionStamperTest, SetReflective) {
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

TYPED_TEST(CFEMDiffusionStamperTest, ConstructorWithReflective) {
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

TYPED_TEST(CFEMDiffusionStamperTest, ConstructorWithReflectiveMap) {
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

/* =============================================================================
 *
 * CFEMDiffusionStamperMPITests
 *
 * Tests using Dealii TestDomains in 1D/2D/3D just to get the cells, matrices,
 * and sparsity patterns.
 *
 * ============================================================================
 */

template <typename TestDomain>
class CFEMDiffusionStamperMPITests
    : public TestDomain,
      public CFEMDiffusionStamperTest<std::integral_constant<int, TestDomain::dimension>> {
 protected:
  static constexpr int dim = TestDomain::dimension;
  dealii::PETScWrappers::MPI::SparseMatrix& system_matrix_ = TestDomain::matrix_1;
  dealii::PETScWrappers::MPI::SparseMatrix& index_hits_    = TestDomain::matrix_2;
  dealii::PETScWrappers::MPI::SparseMatrix& boundary_hits_ = TestDomain::matrix_3;

  dealii::PETScWrappers::MPI::Vector& system_rhs_        = TestDomain::vector_1;
  dealii::PETScWrappers::MPI::Vector& index_hits_vector_ = TestDomain::vector_2;

  void SetUp() override;
};

/* Sets up tests to verify various stamping functions. This includes setting
 * default behavior of the two dependency pointers inherited from
 * CFEMDiffusinStamperTest, and setting up index_hits_ and index_hits_vector_.
 * The values of that matrix and vector, respectively, will indicate how many
 * times that index (or two indices for the matrix) should be stamped if all
 * cells are stamped. This should be identical to if every cell stamps nothing
 * but full matrices or vectors of 1s.
 */
template <typename TestDomain>
void CFEMDiffusionStamperMPITests<TestDomain>::SetUp() {
  auto &mock_definition_ptr = this->mock_definition_ptr;
  // Set up the inherited base class (creates mock dependencies)
  CFEMDiffusionStamperTest<std::integral_constant<int, dim>>::SetUp();
  // Set up the inherited test domain
  TestDomain::SetUpDealii();

  // Set all cell material IDs
  for (const auto& cell : this->cells_) {
    int mat_id = btest::RandomDouble(0, 10);
    cell->set_material_id(mat_id);
  }

  // Set default behaviors
  ON_CALL(*mock_definition_ptr, Cells())
      .WillByDefault(Return(this->cells_));

  dealii::FullMatrix<double> cell_matrix(this->fe_.dofs_per_cell,
                                         this->fe_.dofs_per_cell);
  ON_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillByDefault(Return(cell_matrix));

  dealii::Vector<double> cell_vector(this->fe_.dofs_per_cell);
  ON_CALL(*mock_definition_ptr, GetCellVector())
      .WillByDefault(Return(cell_vector));

  std::vector<dealii::types::global_dof_index> local_dof_indices(this->fe_.dofs_per_cell);

  for (auto cell : this->cells_) {
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      index_hits_vector_(index_i) += 1;
      for (auto index_j : local_dof_indices) {
        index_hits_.add(index_i, index_j, 1);
      }
    }
  }
  index_hits_.compress(dealii::VectorOperation::add);
  index_hits_vector_.compress(dealii::VectorOperation::add);
}

TYPED_TEST_CASE(CFEMDiffusionStamperMPITests, bart::testing::DealiiTestDomains);

TYPED_TEST(CFEMDiffusionStamperMPITests, StampStreaming) {
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;

  int group_number = 1;

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
                FillCellStreamingTerm(_, _, cell, group_number))
        .WillOnce(Invoke(FillMatrixWithOnes));
  }
  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampStreamingTerm(this->system_matrix_, group_number);

  EXPECT_TRUE(CompareMPIMatrices(this->system_matrix_, this->index_hits_));
}

TYPED_TEST(CFEMDiffusionStamperMPITests, StampCollision) {
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;

  int group_number = 1;

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
                FillCellCollisionTerm(_, _, cell, group_number))
        .WillOnce(Invoke(FillMatrixWithOnes));
  }
  EXPECT_CALL(*mock_definition_ptr, GetCellMatrix())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampCollisionTerm(this->system_matrix_, group_number);

  EXPECT_TRUE(CompareMPIMatrices(this->system_matrix_, this->index_hits_));
}

TYPED_TEST(CFEMDiffusionStamperMPITests, StampFixedSource) {
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;
  int group_number = 1;

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
                FillCellFixedSource(_, cell, group_number))
        .WillOnce(Invoke(FillVectorWithOnes3));
  }

  EXPECT_CALL(*mock_definition_ptr, GetCellVector())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampFixedSource(this->system_rhs_, group_number);
  EXPECT_TRUE(CompareMPIVectors(this->index_hits_vector_, this->system_rhs_));
}

TYPED_TEST(CFEMDiffusionStamperMPITests, StampFissionSource) {
  auto& mock_definition_ptr = this->mock_definition_ptr;
  auto& mock_diffusion_ptr = this->mock_diffusion_ptr;
  int group_number = 1;
  double k_effective = 1.04;
  /* Moments for call; it doesn't matter if the moments are uninitialized and
   * empty. We will just be using them as dummys that are caught by the mock, we
   * are just ensuring that the correct things are passed through.
   */
  data::MomentVector in_group_moment;
  data::MomentsMap group_moments;

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*mock_diffusion_ptr,
                FillCellFissionSource(_, cell, group_number, k_effective,
                                      Ref(in_group_moment), Ref(group_moments)))
        .WillOnce(Invoke(FillVectorWithOnes6));
  }

  EXPECT_CALL(*mock_definition_ptr, GetCellVector())
      .WillOnce(DoDefault());

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampFissionSource(this->system_rhs_,
                                  group_number,
                                  k_effective,
                                  in_group_moment,
                                  group_moments);

  EXPECT_TRUE(CompareMPIVectors(this->index_hits_vector_, this->system_rhs_));
}

/* =============================================================================
 *
 * CFEMDiffusionStamperBoundaryMPITests
 *
 * Tests using Dealii TestDomains in 1D/2D/3D just to get the cells, matrices,
 * and sparsity patterns.
 *
 * ============================================================================
 */

template <typename TestDomain>
class CFEMDiffusionStamperBoundaryMPITests :
    public CFEMDiffusionStamperMPITests<TestDomain> {
 protected:
  void SetUp() override;
  void SetUpBoundaries();
};

template <typename TestDomain>
void CFEMDiffusionStamperBoundaryMPITests<TestDomain>::SetUp() {
  CFEMDiffusionStamperMPITests<TestDomain>::SetUp();
  SetUpBoundaries();

  int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
  std::vector<dealii::types::global_dof_index> local_dof_indices(this->fe_.dofs_per_cell);

  for (auto &cell : this->cells_) {
    if (cell->at_boundary()) {
      cell->get_dof_indices(local_dof_indices);
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {
          for (auto index_i : local_dof_indices) {
            for (auto index_j : local_dof_indices) {
              this->boundary_hits_.add(index_i, index_j, 1);
            }
          }
        }
      }
    }
  }
  this->boundary_hits_.compress(dealii::VectorOperation::add);
}

template <typename TestDomain>
void CFEMDiffusionStamperBoundaryMPITests<TestDomain>::SetUpBoundaries() {
  using Boundary = bart::problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
  double zero_tol = 1.0e-14;

  for (auto &cell : this->cells_) {
    for (int face_id = 0; face_id < faces_per_cell; ++face_id) {
      auto face = cell->face(face_id);
      dealii::Point<this->dim> face_center = face->center();

      switch (this->dim) {
        case 3: {
          if (std::fabs(face_center[2]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kZMin));
            break;
          } else if (std::fabs(face_center[2] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kZMax));
            break;
          }
          [[fallthrough]];
        }
        case 2: {
          if (std::fabs(face_center[1]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kYMin));
            break;
          } else if (std::fabs(face_center[1] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kYMax));
            break;
          }
          [[fallthrough]];
        }
        case 1: {
          if (std::fabs(face_center[0]) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kXMin));
            break;
          } else if (std::fabs(face_center[0] - 1) < zero_tol) {
            face->set_boundary_id(static_cast<int>(Boundary::kXMax));
            break;
          }
        }
      }
    }
  }
}

TYPED_TEST_CASE(CFEMDiffusionStamperBoundaryMPITests,
                bart::testing::DealiiTestDomains);

TYPED_TEST(CFEMDiffusionStamperBoundaryMPITests, StampVacuumBoundaryTerm) {
  using BoundaryType =
      typename formulation::scalar::CFEM_DiffusionI<this->dim>::BoundaryType;
  int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
  auto &mock_diffusion_ptr = this->mock_diffusion_ptr;
  auto &mock_definition_ptr = this->mock_definition_ptr;

  for (auto const& cell : this->cells_) {
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

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampBoundaryTerm(this->system_matrix_);

  EXPECT_TRUE(CompareMPIMatrices(this->boundary_hits_, this->system_matrix_));
}

TYPED_TEST(CFEMDiffusionStamperBoundaryMPITests, StampVacuumReflectiveBoundaryTerm) {
  using BoundaryType =
      typename formulation::scalar::CFEM_DiffusionI<this->dim>::BoundaryType;
  using Boundary = problem::Boundary;
  int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
  auto &mock_diffusion_ptr = this->mock_diffusion_ptr;
  auto &mock_definition_ptr = this->mock_definition_ptr;

  for (auto const& cell : this->cells_) {
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

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.AddReflectiveBoundary(Boundary::kXMin);

  test_stamper.StampBoundaryTerm(this->system_matrix_);
}

TYPED_TEST(CFEMDiffusionStamperBoundaryMPITests, StampVacuumAndStreaming) {
  auto &mock_diffusion_ptr = this->mock_diffusion_ptr;
  auto &mock_definition_ptr = this->mock_definition_ptr;

  ON_CALL(*mock_diffusion_ptr, FillCellStreamingTerm(_, _, _, _))
      .WillByDefault(Invoke(FillMatrixWithOnes));
  ON_CALL(*mock_diffusion_ptr, FillBoundaryTerm(_, _, _, _, _))
      .WillByDefault(Invoke(FillMatrixWithOnesBoundary));

  formulation::CFEM_DiffusionStamper<this->dim> test_stamper(
      std::move(mock_diffusion_ptr),
      std::move(mock_definition_ptr));

  test_stamper.StampStreamingTerm(this->system_matrix_, 0);
  test_stamper.StampBoundaryTerm(this->system_matrix_);

  this->index_hits_.add(1, this->boundary_hits_);

  EXPECT_TRUE(CompareMPIMatrices(this->index_hits_, this->system_matrix_));

}

} // namespace