#include "formulation/cfem_saaf_stamper.h"

#include "domain/tests/definition_mock.h"
#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "test_helpers/dealii_test_domain.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::Return, ::testing::_;

template <typename DimensionWrapper>
class CFEM_SAAF_StamperTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using FormulationType =
      formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using DefinitionType = typename domain::DefinitionMock<dim>;
  using InitTokenType = typename formulation::angular::CFEMSelfAdjointAngularFluxI<dim>::InitializationToken;

  std::unique_ptr<FormulationType> formulation_ptr_;
  std::shared_ptr<DefinitionType> definition_ptr_;

  // Mock return objects
  InitTokenType return_token_;

  void SetUp() override;
};
template<typename DimensionWrapper>
void CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp() {
  formulation_ptr_ = std::make_unique<FormulationType>();
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
 */

void StampMatrix(formulation::FullMatrix& to_stamp) {
  for (unsigned int i = 0; i < to_stamp.n_rows(); ++i) {
    for (unsigned int j = 0; j < to_stamp.n_cols(); ++j) {
      to_stamp(i, j) += 1;
    }
  }
}

template <typename DimensionWrapper>
class CFEMSAAFStamperMPITests :
    public CFEM_SAAF_StamperTest<DimensionWrapper>,
    bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  void SetUp() override;

  bart::system::MPISparseMatrix& system_matrix = this->matrix_1;
  bart::system::MPISparseMatrix& index_hits_ = this->matrix_2;
  bart::system::MPIVector& index_hits_vector_ = this->vector_2;
};

template <typename DimensionWrapper>
void CFEMSAAFStamperMPITests<DimensionWrapper>::SetUp() {
  CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp();
  this->SetUpDealii();

  int cell_dofs = this->fe_.dofs_per_cell;

  for (const auto& cell : this->cells_) {
    int material_id = btest::RandomDouble(0, 10);
    cell->set_material_id(material_id);
    std::vector<dealii::types::global_dof_index> local_dof_indices(cell_dofs);
    for (auto index_i : local_dof_indices) {
      index_hits_vector_(index_i) += 1;
      for (auto index_j : local_dof_indices) {
        index_hits_.add(index_i, index_j, 1);
      }
    }
  }
  index_hits_.compress(dealii::VectorOperation::add);
  index_hits_vector_.compress(dealii::VectorOperation::add);

  ON_CALL(*this->definition_ptr_, Cells()).WillByDefault(Return(this->cells_));

  formulation::FullMatrix cell_matrix(cell_dofs, cell_dofs);
  ON_CALL(*this->definition_ptr_, GetCellMatrix())
      .WillByDefault(Return(cell_matrix));
}

TYPED_TEST_SUITE(CFEMSAAFStamperMPITests, bart::testing::AllDimensions);

TYPED_TEST(CFEMSAAFStamperMPITests, StampStreaming) {

}

} // namespace