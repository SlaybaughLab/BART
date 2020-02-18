#include "formulation/cfem_saaf_stamper.h"

#include "domain/tests/definition_mock.h"
#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "quadrature/tests/quadrature_point_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::Return, ::testing::_, ::testing::WithArg,
::testing::Invoke, ::testing::Ref, ::testing::NiceMock;

using bart::testing::CompareMPIMatrices, bart::testing::CompareMPIVectors;

template <typename DimensionWrapper>
class CFEM_SAAF_StamperTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using FormulationType =
      formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using DefinitionType = NiceMock<typename domain::DefinitionMock<dim>>;
  using InitTokenType = typename formulation::angular::CFEMSelfAdjointAngularFluxI<dim>::InitializationToken;

  std::unique_ptr<FormulationType> formulation_ptr_;
  std::shared_ptr<DefinitionType> definition_ptr_;

  // Mock return objects
  InitTokenType return_token_;

  // Observation Pointer
  FormulationType* formulation_obs_ptr_;

  void SetUp() override;
};
template<typename DimensionWrapper>
void CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp() {
  formulation_ptr_ = std::make_unique<FormulationType>();
  formulation_obs_ptr_ = formulation_ptr_.get();
  definition_ptr_ = std::make_shared<DefinitionType>();

  formulation::CellPtr<dim> test_cell_ptr_;
  std::vector<formulation::CellPtr<dim>> cells = {};
  cells.push_back(test_cell_ptr_);

  ON_CALL(*definition_ptr_, Cells()).WillByDefault(Return(cells));
  ON_CALL(*formulation_ptr_, Initialize(_)).WillByDefault(Return(return_token_));
}

TYPED_TEST_SUITE(CFEM_SAAF_StamperTest, bart::testing::AllDimensions);

/* Verify constructor takes and stores dependencies, gets needed variables from
 * dependencies and initializes formulation */
TYPED_TEST(CFEM_SAAF_StamperTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = typename
      formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using DefinitionType = typename domain::DefinitionMock<dim>;

  std::unique_ptr<formulation::CFEM_SAAF_Stamper<dim>> test_stamper;

  EXPECT_CALL(*this->definition_ptr_, Cells()).WillOnce(DoDefault());
  EXPECT_CALL(*this->formulation_ptr_, Initialize(_))
      .WillOnce(DoDefault());
  EXPECT_NO_THROW({
    test_stamper = std::make_unique<formulation::CFEM_SAAF_Stamper<dim>>(
        std::move(this->formulation_ptr_), this->definition_ptr_);
  });

  auto returned_formulation_ptr = test_stamper->formulation_ptr();
  auto returned_definition_ptr = test_stamper->definition_ptr();

  ASSERT_NE(nullptr, returned_formulation_ptr);
  ASSERT_NE(nullptr, returned_definition_ptr);
  ASSERT_NE(nullptr, dynamic_cast<FormulationType*>(returned_formulation_ptr));
  ASSERT_NE(nullptr, dynamic_cast<DefinitionType *>(returned_definition_ptr));
}

/* =============================================================================
 *
 * CFEMSAAFStamperMPITests
 *
 * Tests using deal.ii test domains in 1/2/3D, verifies stamping occurs
 * correctly when a real domain is used.
 *
 * =============================================================================
 */

// HELPER FUNCTIONS ============================================================

void StampMatrix(formulation::FullMatrix& to_stamp) {
  for (unsigned int i = 0; i < to_stamp.n_rows(); ++i) {
    for (unsigned int j = 0; j < to_stamp.n_cols(); ++j) {
      to_stamp(i, j) += 1;
    }
  }
}

void StampVector(formulation::Vector& to_stamp) {
  for (unsigned int i = 0; i < to_stamp.size(); ++i) {
    to_stamp(i) += 1;
  }
}

// TEST SETUP ==================================================================

template <typename DimensionWrapper>
class CFEMSAAFStamperMPITests :
    public CFEM_SAAF_StamperTest<DimensionWrapper>,
    public bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  void SetUp() override;
  static constexpr int dim = DimensionWrapper::value;

  std::unique_ptr<formulation::CFEM_SAAF_Stamper<dim>> test_stamper_;

  // Other test parameters
  std::shared_ptr<quadrature::QuadraturePointI<dim>> quadrature_point_ptr_;
  bart::system::moments::MomentVector in_group_moment_;
  bart::system::moments::MomentsMap moments_map_;

  bart::system::MPISparseMatrix& system_matrix = this->matrix_1;
  bart::system::MPISparseMatrix& index_hits_ = this->matrix_2;
  bart::system::MPISparseMatrix& boundary_hits_ = this->matrix_3;
  bart::system::MPIVector& system_rhs_vector = this->vector_1;
  bart::system::MPIVector& index_hits_vector_ = this->vector_2;
};

template <typename DimensionWrapper>
void CFEMSAAFStamperMPITests<DimensionWrapper>::SetUp() {
  CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp();
  this->SetUpDealii();

  int cell_dofs = this->fe_.dofs_per_cell;

  for (const auto& cell : this->cells_) {
    int material_id = test_helpers::RandomDouble(0, 10);
    cell->set_material_id(material_id);
    std::vector<dealii::types::global_dof_index> local_dof_indices(cell_dofs);
    cell->get_dof_indices(local_dof_indices);
    for (auto index_i : local_dof_indices) {
      index_hits_vector_(index_i) += 1;
      for (auto index_j : local_dof_indices) {
        index_hits_.add(index_i, index_j, 1);
      }
    }
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
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
  index_hits_.compress(dealii::VectorOperation::add);
  index_hits_vector_.compress(dealii::VectorOperation::add);
  boundary_hits_.compress(dealii::VectorOperation::add);

  // Cell matrix and vector
  formulation::FullMatrix cell_matrix(cell_dofs, cell_dofs);
  ON_CALL(*this->definition_ptr_, GetCellMatrix())
      .WillByDefault(Return(cell_matrix));
  formulation::Vector cell_vector(cell_dofs);
  ON_CALL(*this->definition_ptr_, GetCellVector())
      .WillByDefault(Return(cell_vector));

  EXPECT_CALL(*this->formulation_obs_ptr_, Initialize(_)).WillOnce(DoDefault());
  EXPECT_CALL(*this->definition_ptr_, Cells()).WillOnce(Return(this->cells_));

  test_stamper_ = std::make_unique<formulation::CFEM_SAAF_Stamper<dim>>(
      std::move(this->formulation_ptr_), this->definition_ptr_);

  quadrature_point_ptr_ =
      std::make_shared<quadrature::QuadraturePointMock<dim>>();
}

TYPED_TEST_SUITE(CFEMSAAFStamperMPITests, bart::testing::AllDimensions);

// TESTS =======================================================================

TYPED_TEST(CFEMSAAFStamperMPITests, StampBoundary) {
  for (auto const& cell : this->cells_) {
    if (cell->at_boundary()) {
      int faces_per_cell = dealii::GeometryInfo<this->dim>::faces_per_cell;
      for (int face = 0; face < faces_per_cell; ++face) {
        if (cell->face(face)->at_boundary()) {

          EXPECT_CALL(*this->formulation_obs_ptr_,
                      FillBoundaryBilinearTerm(_,_,cell, domain::FaceIndex(face),
                                               this->quadrature_point_ptr_,
                                               system::EnergyGroup(1)))
              .WillOnce(::testing::WithArg<0>(::testing::Invoke(StampMatrix)));
        }
      }
    }
  }

  EXPECT_CALL(*this->definition_ptr_, GetCellMatrix()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
                    this->test_stamper_->StampBoundaryBilinearTerm(
                        this->system_matrix,
                        this->quadrature_point_ptr_,
                        system::EnergyGroup(1));
                  });

  EXPECT_TRUE(CompareMPIMatrices(this->boundary_hits_, this->system_matrix));

}

TYPED_TEST(CFEMSAAFStamperMPITests, StampCollision) {

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
        FillCellCollisionTerm(_,_,cell,system::EnergyGroup(1)))
        .WillOnce(::testing::WithArg<0>(::testing::Invoke(StampMatrix)));
  }

  EXPECT_CALL(*this->definition_ptr_, GetCellMatrix()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
                    this->test_stamper_->StampCollisionTerm(
                        this->system_matrix, system::EnergyGroup(1));
                  });

  EXPECT_TRUE(CompareMPIMatrices(this->index_hits_, this->system_matrix));
}

TYPED_TEST(CFEMSAAFStamperMPITests, StampFissionSource) {
  const double k_effective = 1.14;

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
        FillCellFissionSourceTerm(_,_,cell, this->quadrature_point_ptr_,
                                  system::EnergyGroup(1), k_effective,
                                  Ref(this->in_group_moment_),
                                  Ref(this->moments_map_)))
        .WillOnce(WithArg<0>(Invoke(StampVector)));
  }

  EXPECT_CALL(*this->definition_ptr_, GetCellVector()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
    this->test_stamper_->StampFissionSourceTerm(
        this->system_rhs_vector,
        this->quadrature_point_ptr_,
        system::EnergyGroup(1),
        k_effective,
        this->in_group_moment_,
        this->moments_map_);
                  });

  EXPECT_TRUE(CompareMPIVectors(this->index_hits_vector_,
                                this->system_rhs_vector));

}

TYPED_TEST(CFEMSAAFStamperMPITests, StampFixedSource) {
  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellFixedSourceTerm(_, _, cell, this->quadrature_point_ptr_,
                    system::EnergyGroup(1)))
        .WillOnce(WithArg<0>(Invoke(StampVector)));
  }
  EXPECT_CALL(*this->definition_ptr_, GetCellVector()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
    this->test_stamper_->StampFixedSourceTerm(
        this->system_rhs_vector,
        this->quadrature_point_ptr_,
        system::EnergyGroup(1));
  });

  EXPECT_TRUE(CompareMPIVectors(this->system_rhs_vector,
                                this->index_hits_vector_));
}

TYPED_TEST(CFEMSAAFStamperMPITests, StampScatteringSource) {

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellScatteringSourceTerm(_,_,cell,
                                             this->quadrature_point_ptr_,
                                             system::EnergyGroup(1),
                                             Ref(this->in_group_moment_),
                                             Ref(this->moments_map_)))
        .WillOnce(WithArg<0>(Invoke(StampVector)));
  }

  EXPECT_CALL(*this->definition_ptr_, GetCellVector()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
                    this->test_stamper_->StampScatteringSourceTerm(
                        this->system_rhs_vector,
                        this->quadrature_point_ptr_,
                        system::EnergyGroup(1),
                        this->in_group_moment_,
                        this->moments_map_);
                  });

  EXPECT_TRUE(CompareMPIVectors(this->index_hits_vector_,
                                this->system_rhs_vector));

}

TYPED_TEST(CFEMSAAFStamperMPITests, StampStreaming) {

  for (auto const& cell : this->cells_) {
    EXPECT_CALL(*this->formulation_obs_ptr_,
                FillCellStreamingTerm(_,_,cell,
                                      this->quadrature_point_ptr_,
                                      system::EnergyGroup(1)))
        .WillOnce(::testing::WithArg<0>(::testing::Invoke(StampMatrix)));
  }

  EXPECT_CALL(*this->definition_ptr_, GetCellMatrix()).WillOnce(DoDefault());

  EXPECT_NO_THROW({
                    this->test_stamper_->StampStreamingTerm(
                        this->system_matrix,
                        this->quadrature_point_ptr_,
                        system::EnergyGroup(1));
                  });

  EXPECT_TRUE(CompareMPIMatrices(this->index_hits_, this->system_matrix));
}

} // namespace