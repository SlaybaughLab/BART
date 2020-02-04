#include "iteration/updater/source_updater_gauss_seidel.h"

#include <memory>
#include <iteration/updater/angular_source_updater_gauss_seidel.h>

#include "test_helpers/gmock_wrapper.h"
#include "quadrature/tests/quadrature_set_mock.h"
#include "formulation/tests/angular_stamper_mock.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class IterationUpdaterAngularSourceUpdaterGaussSeidelTest :
    public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;

  std::shared_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  std::unique_ptr<UpdaterType> test_updater_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationUpdaterAngularSourceUpdaterGaussSeidelTest<DimensionWrapper>::SetUp() {
  stamper_ptr_ = std::make_shared<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  test_updater_ = std::make_unique<UpdaterType>(stamper_ptr_, quadrature_set_ptr_);
}

TYPED_TEST_SUITE(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest, Constructor) {
  constexpr int dim = this->dim;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  EXPECT_NO_THROW({ UpdaterType test_updater(stamper_ptr, quadrature_set_ptr); });
  EXPECT_ANY_THROW({ UpdaterType test_updater(nullptr, quadrature_set_ptr); });
  EXPECT_ANY_THROW({ UpdaterType test_updater(stamper_ptr, nullptr); });
}

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest, Getters) {
  auto quadrature_set_ptr = this->test_updater_->quadrature_set_ptr();
  auto stamper_ptr = this->test_updater_->stamper_ptr();

  EXPECT_EQ(quadrature_set_ptr, this->quadrature_set_ptr_.get());
  EXPECT_EQ(stamper_ptr, this->stamper_ptr_.get());
}

} // namespace
