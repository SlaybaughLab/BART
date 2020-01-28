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
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using FixedUpdaterType = iteration::updater::AngularFixedUpdater<formulation::AngularStamperI<dim>>;

  std::shared_ptr<StamperType> stamper_ptr_;
  std::shared_ptr<QuadratureSetType> quadrature_set_ptr_;
  std::unique_ptr<FixedUpdaterType> test_updater_ptr_;

  void SetUp() override;
};

template <typename DimensionWrapper>
void IterationUpdaterAngularFixedUpdaterTest<DimensionWrapper>::SetUp() {
  stamper_ptr_ = std::make_shared<StamperType>();
  quadrature_set_ptr_ = std::make_shared<QuadratureSetType>();

  test_updater_ptr_ = std::make_unique<FixedUpdaterType>(
      stamper_ptr_, quadrature_set_ptr_);
}

TYPED_TEST_SUITE(IterationUpdaterAngularFixedUpdaterTest,
                 bart::testing::AllDimensions);

TYPED_TEST(IterationUpdaterAngularFixedUpdaterTest, Constructor) {
  constexpr int dim = this->dim;
  using StamperType = formulation::AngularStamperMock<dim>;
  using QuadratureSetType = quadrature::QuadratureSetMock<dim>;
  using FixedUpdaterType =
      iteration::updater::AngularFixedUpdater<formulation::AngularStamperI<dim>>;

  auto stamper_ptr = std::make_shared<StamperType>();
  auto quadrature_set_ptr = std::make_shared<QuadratureSetType>();

  EXPECT_NO_THROW({
    FixedUpdaterType test_updater(stamper_ptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, quadrature_set_ptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, nullptr);
  });

  EXPECT_ANY_THROW({
    FixedUpdaterType test_updater(nullptr, quadrature_set_ptr);
  });
}

TYPED_TEST(IterationUpdaterAngularFixedUpdaterTest, Getters) {
  auto stamper_ptr = this->test_updater_ptr_->stamper_ptr();
  auto quadrature_set_ptr = this->test_updater_ptr_->quadrature_set_ptr();

  ASSERT_NE(nullptr, stamper_ptr);
  EXPECT_EQ(stamper_ptr, this->stamper_ptr_.get());

  ASSERT_NE(nullptr, quadrature_set_ptr);
  EXPECT_EQ(quadrature_set_ptr, this->quadrature_set_ptr_.get());
}



} // namespace
