#include "formulation/cfem_saaf_stamper.h"

#include "domain/tests/definition_mock.h"
#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class CFEM_SAAF_StamperTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using FormulationType =
      formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using DefinitionType = typename domain::DefinitionMock<dim>;

  std::unique_ptr<FormulationType> formulation_ptr_;
  std::shared_ptr<DefinitionType> definition_ptr_;
  void SetUp() override;
};
template<typename DimensionWrapper>
void CFEM_SAAF_StamperTest<DimensionWrapper>::SetUp() {
  formulation_ptr_ = std::make_unique<FormulationType>();
  definition_ptr_ = std::make_shared<DefinitionType>();
}

TYPED_TEST_SUITE(CFEM_SAAF_StamperTest, bart::testing::AllDimensions);

TYPED_TEST(CFEM_SAAF_StamperTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = typename
      formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using DefinitionType = typename domain::DefinitionMock<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  std::unique_ptr<formulation::CFEM_SAAF_Stamper<dim>> test_stamper;

  EXPECT_NO_THROW({
    test_stamper = std::make_unique<formulation::CFEM_SAAF_Stamper<dim>>(
        std::move(formulation_ptr), this->definition_ptr_);
  });

  auto returned_formulation_ptr = test_stamper->formulation_ptr();
  auto returned_definition_ptr = test_stamper->definition_ptr();

  ASSERT_NE(nullptr, returned_formulation_ptr);
  ASSERT_NE(nullptr, returned_definition_ptr);
  ASSERT_NE(nullptr, dynamic_cast<FormulationType*>(returned_formulation_ptr));
  ASSERT_NE(nullptr, dynamic_cast<DefinitionType *>(returned_definition_ptr));
}

} // namespace