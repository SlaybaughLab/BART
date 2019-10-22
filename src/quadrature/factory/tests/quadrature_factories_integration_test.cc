#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/angular/scalar_angular.h"
#include "quadrature/angular/level_symmetric_gaussian.h"
#include "quadrature/ordinate.h"
#include "quadrature/quadrature_point.h"

#include "test_helpers/gmock_wrapper.h"


namespace  {

using namespace bart;

/* These tests verify the operation of the functions in quadrature::factory
 * that are designed to instantiate the different implementations of all the
 * classes in the quadrature namespace. As they actually instantiate classes,
 * they are integration tests, as they require the instantiated classes to be
 * implemented in some form.
 */
template <typename DimensionWrapper>
class QuadratureFactoriesIntegrationTest : public ::testing::Test {
 public:
  static constexpr int dim = DimensionWrapper::value;
};

TYPED_TEST_CASE(QuadratureFactoriesIntegrationTest, bart::testing::AllDimensions);

/* Call to MakeOrdinate specifying default implementation of OrdinateI should
 * return the correct type.
 */
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeOrdinate) {
  constexpr int dim = this->dim;
  auto ordinate_ptr = quadrature::factory::MakeOrdinatePtr<dim>(
      quadrature::OrdinateType::kDefault);

  ASSERT_NE(nullptr, ordinate_ptr);

  using ExpectedType = quadrature::Ordinate<dim>;
  EXPECT_NE(nullptr, dynamic_cast<ExpectedType*>(ordinate_ptr.get()));
}

/* Call to MakeQuadraturePoint specifying default implementation of
 * QuadraturePoint should return the correct type.
 */
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeQuadraturePoint) {
  constexpr int dim = this->dim;
  auto quadrature_point_ptr = quadrature::factory::MakeQuadraturePointPtr<dim>(
      quadrature::QuadraturePointImpl::kDefault);

  ASSERT_NE(nullptr, quadrature_point_ptr);

  using ExpectedType = quadrature::QuadraturePoint<dim>;
  EXPECT_NE(nullptr, dynamic_cast<ExpectedType*>(quadrature_point_ptr.get()));
}

/* Call to MakeAngularQuadratureGeneratorPtr specifying no type should return
 * a runtime error. */
TYPED_TEST(QuadratureFactoriesIntegrationTest,
    MakeAngularQuadratureGenTestNoType) {
  constexpr int dim = this->dim;
  const int order_value = 4;

  EXPECT_ANY_THROW({
    quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
        quadrature::Order(order_value),
        quadrature::AngularQuadratureSetType::kNone);
  });
}

/* Call to MakeAngularQuadratureGeneratorPtr specifying a scalar quadrature
 * should return the correct type. */
TYPED_TEST(QuadratureFactoriesIntegrationTest,
           MakeAngularQuadratureGenTestScalar) {
  constexpr int dim = this->dim;
  const int order_value = 4;

  auto quadrature_generator_ptr =
      quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
          quadrature::Order(order_value),
          quadrature::AngularQuadratureSetType::kScalar);

  ASSERT_NE(nullptr, quadrature_generator_ptr);
  using ExpectedType = quadrature::angular::ScalarAngular<dim>;
  EXPECT_NE(nullptr,
      dynamic_cast<ExpectedType*>(quadrature_generator_ptr.get()));
  EXPECT_EQ(quadrature_generator_ptr->order(), 0);
}

/* Call to MakeAngularQuadratureGeneratorPtr specifying a level symmetric gaussian
 * should return the correct type. If dim !=3, should throw an error */
TYPED_TEST(QuadratureFactoriesIntegrationTest,
           MakeAngularQuadratureGenTestLSGaussian) {
  constexpr int dim = this->dim;
  const int order_value = 4;
  if (dim == 3) {
        auto quadrature_generator_ptr =
        quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
            quadrature::Order(order_value),
            quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);

    ASSERT_NE(nullptr, quadrature_generator_ptr);
    using ExpectedType = quadrature::angular::LevelSymmetricGaussian;
    EXPECT_NE(nullptr,
              dynamic_cast<ExpectedType*>(quadrature_generator_ptr.get()));
    EXPECT_EQ(quadrature_generator_ptr->order(), order_value);
  } else {
    EXPECT_ANY_THROW({
      quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
          quadrature::Order(order_value),
          quadrature::AngularQuadratureSetType::kLevelSymmetricGaussian);
    });
  }
}

} // namespace