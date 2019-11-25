#include "formulation/angular/cfem_self_adjoint_angular_flux.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/dealii_test_domain.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class FormulationAngularCFEMSelfAdjointAngularFluxTest :
    public ::testing::Test ,
    bart::testing::DealiiTestDomain<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationAngularCFEMSelfAdjointAngularFluxTest<DimensionWrapper>::SetUp() {
  this->SetUpDealii();
}

TYPED_TEST_CASE(FormulationAngularCFEMSelfAdjointAngularFluxTest,
                bart::testing::AllDimensions);

TYPED_TEST(FormulationAngularCFEMSelfAdjointAngularFluxTest, Dummy) {
  EXPECT_TRUE(true);
}


} // namespace