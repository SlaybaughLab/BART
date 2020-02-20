#include "formulation/updater/saaf_updater.h"

#include "formulation/angular/tests/cfem_self_adjoint_angular_flux_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationUpdaterSAAFTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(FormulationUpdaterSAAFTest, bart::testing::AllDimensions);

TYPED_TEST(FormulationUpdaterSAAFTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;

  EXPECT_NO_THROW({
    test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                     std::move(stamper_ptr));
  });
  EXPECT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->stamper_ptr(), nullptr);
}

TYPED_TEST(FormulationUpdaterSAAFTest, ConstructorBadDepdendencies) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::CFEMSelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;

  for (bool formulation_good : {true, false}) {
    for (bool stamper_good : {true, false}) {
      auto formulation_ptr = formulation_good ?
                             std::make_unique<FormulationType>() : nullptr;
      auto stamper_ptr = stamper_good ?
                         std::make_unique<StamperType>() : nullptr;
      if (!formulation_good || !stamper_good) {
        std::unique_ptr<UpdaterType> test_updater_ptr;
        EXPECT_ANY_THROW({
          test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                           std::move(stamper_ptr));
                        });
      }
    }
  }
}


} // namespace
