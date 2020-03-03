#include "formulation/updater/diffusion_updater.h"

#include "formulation/scalar/tests/diffusion_mock.h"
#include "formulation/tests/stamper_mock.h"
#include "formulation/updater/tests/updater_tests.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_assertions.h"

namespace {

using namespace bart;

template <typename DimensionWrapper>
class FormulationUpdaterDiffusionTest :
    public bart::formulation::updater::test_helpers::UpdaterTests<DimensionWrapper::value> {
 public:
  static constexpr int dim = DimensionWrapper::value;

  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  std::unique_ptr<UpdaterType> test_updater_ptr_;

  FormulationType* formulation_obs_ptr_;
  StamperType* stamper_obs_ptr_;

  void SetUp() override;
};

TYPED_TEST_SUITE(FormulationUpdaterDiffusionTest, bart::testing::AllDimensions);

template <typename DimensionWrapper>
void FormulationUpdaterDiffusionTest<DimensionWrapper>::SetUp() {
  bart::formulation::updater::test_helpers::UpdaterTests<dim>::SetUp();
  auto formulation_ptr = std::make_unique<FormulationType>();
  formulation_obs_ptr_ = formulation_ptr.get();
  auto stamper_ptr = this->MakeStamper();
  stamper_obs_ptr_ = stamper_ptr.get();
}

// ===== CONSTRUCTOR TESTS =====================================================

TYPED_TEST(FormulationUpdaterDiffusionTest, Constructor) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_NO_THROW(
      {
        test_updater_ptr = std::make_unique<UpdaterType>(
            std::move(formulation_ptr), std::move(stamper_ptr));
      }
  );
  ASSERT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  ASSERT_NE(test_updater_ptr->stamper_ptr(), nullptr);
}

TYPED_TEST(FormulationUpdaterDiffusionTest, ConstructorReflective) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unordered_set<problem::Boundary> reflective_boundaries{problem::Boundary::kXMin};
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_NO_THROW(
      {
        test_updater_ptr = std::make_unique<UpdaterType>(
            std::move(formulation_ptr), std::move(stamper_ptr),
            reflective_boundaries);
      }
  );
  ASSERT_NE(test_updater_ptr->formulation_ptr(), nullptr);
  ASSERT_NE(test_updater_ptr->stamper_ptr(), nullptr);
  EXPECT_EQ(reflective_boundaries, test_updater_ptr->reflective_boundaries());
}

TYPED_TEST(FormulationUpdaterDiffusionTest, ConstructorBadDependencies) {
  constexpr int dim = this->dim;
  using FormulationType = formulation::scalar::DiffusionMock<dim>;
  using StamperType = formulation::StamperMock<dim>;
  using UpdaterType = formulation::updater::DiffusionUpdater<dim>;

  auto formulation_ptr = std::make_unique<FormulationType>();
  auto stamper_ptr = std::make_unique<StamperType>();
  std::unique_ptr<UpdaterType> test_updater_ptr;
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
      nullptr, std::move(stamper_ptr)); });
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
        std::move(formulation_ptr), nullptr); });
  EXPECT_ANY_THROW({test_updater_ptr = std::make_unique<UpdaterType>(
      nullptr, nullptr); });
}

} // namespace