#include "iteration/updater/angular_fixed_updater.h"
#include "formulation/tests/angular_stamper_mock.h"
#include "quadrature/tests/quadrature_set_mock.h"

#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

template <typename DimensionWrapper>
class IterationUpdaterAngularFixedUpdaterTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_SUITE(IterationUpdaterAngularFixedUpdaterTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularFixedUpdaterTest, Constructor) {
  constexpr int dim = this->dim;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using FixedUpdaterType = iteration::updater::AngularFixedUpdater<StamperType>;

  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();

  EXPECT_NO_THROW({
    FixedUpdaterType test_updater(stamper_ptr, quadrature_set_ptr);
  });
}

} // namespace
