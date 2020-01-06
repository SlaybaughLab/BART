#include "formulation/cfem_saaf_stamper.h"

#include "domain/tests/definition_mock.h"
#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "test_helpers/gmock_wrapper.h"

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
  std::vector<formulation::CellPtr<dim>> cells_ = {};
  InitTokenType return_token_;

  void SetUp() override;
};
template<typename DimensionWrapper>
void CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp() {
  formulation_ptr_ = std::make_unique<FormulationType>();
  definition_ptr_ = std::make_shared<DefinitionType>();

  formulation::CellPtr<dim> test_cell_ptr_;
  cells_.push_back(test_cell_ptr_);

  ON_CALL(*definition_ptr_, Cells()).WillByDefault(Return(cells_));
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
  EXPECT_CALL(*this->formulation_ptr_, Initialize(this->cells_.at(0)))
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

} // namespace