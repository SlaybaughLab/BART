#include "quadrature/factory/quadrature_factories.h"

#include "quadrature/angular/gauss_legendre.h"
#include "quadrature/angular/level_symmetric_gaussian.h"
#include "quadrature/ordinate.hpp"
#include "quadrature/quadrature_point.hpp"
#include "quadrature/quadrature_set.h"
#include "quadrature/calculators/scalar_moment.h"
#include "quadrature/calculators/spherical_harmonic_zeroth_moment.h"
#include "quadrature/utility/quadrature_utilities.h"

#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"


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

/* Call to MakeAngularQuadratureGeneratorPtr specifying a GaussLegendre
 * should return the correct type. If dim !=1, should throw an error */
TYPED_TEST(QuadratureFactoriesIntegrationTest,
           MakeAngularQuadratureGenTestGaussLegendre) {
  constexpr int dim = this->dim;
  const int order_value = 4;
  if (dim == 1) {
    auto quadrature_generator_ptr =
        quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
            quadrature::Order(order_value),
            quadrature::AngularQuadratureSetType::kGaussLegendre);

    ASSERT_NE(nullptr, quadrature_generator_ptr);
    using ExpectedType = quadrature::angular::GaussLegendre;
    EXPECT_NE(nullptr,
              dynamic_cast<ExpectedType*>(quadrature_generator_ptr.get()));
    EXPECT_EQ(quadrature_generator_ptr->order(), order_value);
  } else {
    EXPECT_ANY_THROW({
                       quadrature::factory::MakeAngularQuadratureGeneratorPtr<dim>(
                           quadrature::Order(order_value),
                           quadrature::AngularQuadratureSetType::kGaussLegendre);
                     });
  }
}



/* Call to MakeQuadratureSet specifying default implementation should return
 * the correct type. */
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeQuadratureSetTest) {
  constexpr int dim = this->dim;

  auto quadrature_set_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>(
      quadrature::QuadratureSetImpl::kDefault);

  ASSERT_NE(nullptr, quadrature_set_ptr);

  using ExpectedType = quadrature::QuadratureSet<dim>;

  ASSERT_NE(nullptr,
            dynamic_cast<ExpectedType*>(quadrature_set_ptr.get()));
}

TYPED_TEST(QuadratureFactoriesIntegrationTest, FillQuadratureSet) {
  const int dim = this->dim;
  const int n_points = 3;
  const int n_quadrants = std::pow(2, dim);

  // Generate random quadrature points
  std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight>>
      quadrature_points;

  for (int i = 0; i < n_points; ++i) {
    auto random_position = test_helpers::RandomVector(dim, 1, 10);
    auto random_weight = test_helpers::RandomDouble(0, 2);
    std::array<double, dim> position;
    for (int j = 0; j < dim; ++j)
      position.at(j) = random_position.at(j);
    quadrature_points.emplace_back(quadrature::CartesianPosition<dim>(position),
                                quadrature::Weight(random_weight));
  }

  // Distribute in all positive X quadrants
  auto distributed_points =
      quadrature::utility::GenerateAllPositiveX<dim>(quadrature_points);

  auto quadrature_set_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

  quadrature::factory::FillQuadratureSet<dim>(quadrature_set_ptr.get(),
                                              distributed_points);

  EXPECT_EQ(quadrature_set_ptr->size(), n_points*n_quadrants);

  for (const auto& quadrature_point_ptr : *quadrature_set_ptr) {

    // Verify this point is in the original set of points if it is all positive,
    // need to construct a pair to do this
    quadrature::CartesianPosition<dim> position(quadrature_point_ptr->ordinate()->cartesian_position());
    quadrature::Weight weight(quadrature_point_ptr->weight());

    if (std::all_of(position.get().cbegin(),
                    position.get().cend(), [](int i) { return i >=0;})) {
      auto point_pair = std::make_pair(position, weight);
      EXPECT_EQ(1, std::count(distributed_points.cbegin(), distributed_points.cend(),
                              point_pair));
    }

    // This point should have a reflection, whose reflection should be this point
    // and the reflection should have the opposite coordinate
    auto reflection_ptr = quadrature_set_ptr->GetReflection(quadrature_point_ptr);
    ASSERT_NE(nullptr, reflection_ptr);
    auto re_reflection_ptr = quadrature_set_ptr->GetReflection(reflection_ptr);
    EXPECT_EQ(re_reflection_ptr, quadrature_point_ptr);
    EXPECT_EQ(quadrature_point_ptr->weight(), reflection_ptr->weight());

    for (int i = 0; i < dim; ++i) {
      EXPECT_EQ(quadrature_point_ptr->ordinate()->cartesian_position().at(i),
                -reflection_ptr->ordinate()->cartesian_position().at(i));
    }
  }
  // Subsequent calls should throw an error (quadrature set is not empty)
  EXPECT_ANY_THROW({
    quadrature::factory::FillQuadratureSet<dim>(quadrature_set_ptr.get(),
                                                distributed_points);
                   });
}

TYPED_TEST(QuadratureFactoriesIntegrationTest, FillQuadratureSetScalar) {
  const int dim = this->dim;
  std::vector<std::pair<quadrature::CartesianPosition<dim>, quadrature::Weight >>
      quadrature_points;

  std::array<double, dim> zero_position;
  zero_position.fill(0);

  quadrature_points.emplace_back(quadrature::CartesianPosition<dim>(zero_position),
                                 quadrature::Weight(1.0));

  auto quadrature_set_ptr = quadrature::factory::MakeQuadratureSetPtr<dim>();

  quadrature::factory::FillQuadratureSet<dim>(quadrature_set_ptr.get(),
                                              quadrature_points);

  ASSERT_EQ(quadrature_set_ptr->size(), 1);
  auto quadrature_point_ptr = *quadrature_set_ptr->begin();

  for (int i = 0; i < dim; ++i) {
    EXPECT_EQ(0, quadrature_point_ptr->ordinate()->cartesian_position().at(i));
  }
}

// MakeMomentCalculator should return the correct scalar moment implementation
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeMomentCalculatorScalar) {
  const int dim = this->dim;

  auto moment_calculator_ptr = quadrature::factory::MakeMomentCalculator<dim>(
      quadrature::MomentCalculatorImpl::kScalarMoment);

  ASSERT_NE(moment_calculator_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<quadrature::calculators::ScalarMoment*>(
                moment_calculator_ptr.get()));
}

// MakeMomentCalculator should return the correct zeroth-only moment implementation
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeMomentCalculatorZeroth) {
  const int dim = this->dim;
  auto quadrature_set_mock = quadrature::factory::MakeQuadratureSetPtr<dim>();

  auto moment_calculator_ptr = quadrature::factory::MakeMomentCalculator<dim>(
      quadrature::MomentCalculatorImpl::kZerothMomentOnly,
      quadrature_set_mock);

  ASSERT_NE(moment_calculator_ptr, nullptr);
  EXPECT_NE(nullptr,
            dynamic_cast<quadrature::calculators::SphericalHarmonicZerothMoment<dim>*>(
                moment_calculator_ptr.get()));
}

// MakeMomentCalculator should throw if quadrature set is null and an
// implementation is requested that requires a quadrature set
TYPED_TEST(QuadratureFactoriesIntegrationTest, MakeMomentCalculatorZerothBadSet) {
  const int dim = this->dim;

  EXPECT_ANY_THROW({
    quadrature::factory::MakeMomentCalculator<dim>(
        quadrature::MomentCalculatorImpl::kZerothMomentOnly, nullptr);
                   });
}




} // namespace