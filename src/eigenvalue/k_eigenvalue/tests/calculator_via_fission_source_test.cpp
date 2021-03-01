#include "eigenvalue/k_eigenvalue/calculator_via_fission_source.hpp"

#include "calculator/cell/tests/total_aggregated_fission_source_mock.hpp"
#include "system/system.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return;

class EigenvalueKEigenvalueCalculatorViaFissionSourceTest : public ::testing::Test {
 public:

  using TotalAggregatedFissionSourceType = calculator::cell::TotalAggregatedFissionSourceMock;
  using MomentsType = system::moments::SphericalHarmonicMock;
  using TestCalculator = eigenvalue::k_eigenvalue::CalculatorViaFissionSource;

  std::unique_ptr<TestCalculator> test_k_calculator_;

  // Supporting objects
  TotalAggregatedFissionSourceType* fission_source_mock_obs_ptr_;
  MomentsType* current_moments_obs_ptr_;
  system::System test_system;

  // test parameters
  const double initial_k_eigenvalue_{ 1.0 };
  const double initial_fission_source_{ 5.0 };

  void SetUp() override;
};

void EigenvalueKEigenvalueCalculatorViaFissionSourceTest::SetUp() {
  auto fission_source_mock_ptr = std::make_unique<TotalAggregatedFissionSourceType>();
  auto current_moments_mock_ptr = std::make_unique<MomentsType>();

  fission_source_mock_obs_ptr_ = fission_source_mock_ptr.get();
  current_moments_obs_ptr_ = current_moments_mock_ptr.get();

  test_k_calculator_ = std::make_unique<TestCalculator>(std::move(fission_source_mock_ptr),
                                                      initial_k_eigenvalue_, initial_fission_source_);
  test_system.current_moments = std::move(current_moments_mock_ptr);
}

TEST_F(EigenvalueKEigenvalueCalculatorViaFissionSourceTest, Constructor) {
  EXPECT_EQ(test_k_calculator_->k_eigenvalue(), std::nullopt);
  EXPECT_EQ(test_k_calculator_->initial_k_eigenvalue(), initial_k_eigenvalue_);
  EXPECT_EQ(test_k_calculator_->initial_fission_source(), initial_fission_source_);
  EXPECT_EQ(test_k_calculator_->current_fission_source(), std::nullopt);
  EXPECT_EQ(test_k_calculator_->fission_source_calculator(), fission_source_mock_obs_ptr_);
}

TEST_F(EigenvalueKEigenvalueCalculatorViaFissionSourceTest, BadConstructor) {

  std::array<double,2> bad_keff{0, -1.3};
  std::array<double,2> bad_fission_sources{0, -10.0};

  for (auto const keff : bad_keff) {
    for (auto const fission_source : bad_fission_sources) {
      EXPECT_ANY_THROW({
        TestCalculator(std::make_unique<TotalAggregatedFissionSourceType>(), keff, fission_source);
      });
    }
  }
}

TEST_F(EigenvalueKEigenvalueCalculatorViaFissionSourceTest, BadFissionSources) {
  EXPECT_CALL(*fission_source_mock_obs_ptr_, AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(-1));
  EXPECT_ANY_THROW(test_k_calculator_->CalculateK_Eigenvalue(test_system));

  EXPECT_CALL(*fission_source_mock_obs_ptr_, AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(0));
  EXPECT_ANY_THROW(test_k_calculator_->CalculateK_Eigenvalue(test_system));
}

TEST_F(EigenvalueKEigenvalueCalculatorViaFissionSourceTest, CalculateKEff) {

  const std::array<double, 3> fission_sources{5.3, 2.6, 1.3};
  const std::array<double, 3> expected_result{ 1.06, 0.52, 0.26 };
  const int steps{ 3 };

  for (int step = 0; step < steps; ++step) {
    EXPECT_CALL(*fission_source_mock_obs_ptr_, AggregatedFissionSource(test_system.current_moments.get()))
        .WillOnce(Return(fission_sources.at(step)));
    const auto k_effective = test_k_calculator_->CalculateK_Eigenvalue(test_system);
    test_system.k_effective = k_effective;

    EXPECT_NEAR(k_effective, expected_result.at(step), 1e-12);
    EXPECT_DOUBLE_EQ(test_k_calculator_->current_fission_source().value_or(0), fission_sources.at(step));
    EXPECT_DOUBLE_EQ(test_k_calculator_->k_eigenvalue().value_or(0), expected_result.at(step));
  }
}

} // namespace