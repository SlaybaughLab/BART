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
};

TYPED_TEST_SUITE(IterationUpdaterAngularSourceUpdaterGaussSeidelTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularSourceUpdaterGaussSeidelTest, Constructor) {
  constexpr int dim = this->dim;
  using UpdaterType = iteration::updater::AngularSourceUpdaterGaussSeidel<formulation::AngularStamperI<dim>>;
  using StamperType = formulation::AngularStamperMock<dim>;
  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<quadrature::QuadratureSetMock<dim>>();

  EXPECT_NO_THROW({
    UpdaterType test_updater(stamper_ptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    UpdaterType test_updater(nullptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    UpdaterType test_updater(stamper_ptr, nullptr);
  });
}

} // namespace
