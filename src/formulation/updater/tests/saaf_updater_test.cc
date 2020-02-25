#include "formulation/updater/saaf_updater.h"

#include "quadrature/tests/quadrature_set_mock.h"
#include "formulation/angular/tests/self_adjoint_angular_flux_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationUpdaterSAAFTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  // Test object
  std::unique_ptr<UpdaterType> test_updater_ptr;

  // Pointers to mocks
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  FormulationType* formulation_obs_ptr_;
  StamperType* stamper_obs_ptr_;
  void SetUp() override;
};

template <typename DimensionWrapper>
void FormulationUpdaterSAAFTest<DimensionWrapper>::SetUp() {
  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();
  test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                   std::move(stamper_ptr),
                                                   quadrature_set_ptr_);
}

TYPED_TEST_SUITE(FormulationUpdaterSAAFTest, bart::testing::AllDimensions);

// ====== Constructor tests ====================================================

TYPED_TEST(FormulationUpdaterSAAFTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;

  EXPECT_NO_THROW({
    test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                                                     std::move(stamper_ptr),
                                                     quadrature_set_ptr);
  });
  EXPECT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->stamper_ptr(), nullptr);
  EXPECT_NE(test_updater_ptr->quadrature_set_ptr(), nullptr);
}

TYPED_TEST(FormulationUpdaterSAAFTest, ConstructorBadDepdendencies) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::angular::SelfAdjointAngularFluxMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::SAAFUpdater<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  for (bool formulation_good : {true, false}) {
    for (bool stamper_good : {true, false}) {
      for (bool quadrature_set_good : {true, false}) {
        auto formulation_ptr = formulation_good ?
                               std::make_unique<FormulationType>() : nullptr;
        auto stamper_ptr = stamper_good ?
                           std::make_unique<StamperType>() : nullptr;
        auto quadrature_set_ptr = quadrature_set_good ?
            std::make_shared<QuadratureSetType>() : nullptr;
        if (!formulation_good || !stamper_good || !quadrature_set_good) {
          std::unique_ptr<UpdaterType> test_updater_ptr;
          EXPECT_ANY_THROW({
            test_updater_ptr = std::make_unique<UpdaterType>(std::move(formulation_ptr),
                std::move(stamper_ptr), quadrature_set_ptr);
          });
        }
      }
    }
  }
}

// ===== Update Fixed Terms Tests
TYPED_TEST(FormulationUpdaterSAAFTest, UpdateFixedTermsTest) {

}

} // namespace
