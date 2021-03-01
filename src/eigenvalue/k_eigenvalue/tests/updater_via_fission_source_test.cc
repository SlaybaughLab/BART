#include <calculator/cell/tests/integrated_fission_source_mock.hpp>
#include "eigenvalue/k_eigenvalue/updater_via_fission_source.h"

#include "calculator/cell/tests/total_aggregated_fission_source_mock.hpp"
#include "system/system.hpp"
#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"

namespace  {

using namespace bart;

using ::testing::Return;

class EigKEffUpdaterViaFissionSourceTest : public ::testing::Test {
 protected:

  using TotalAggregatedFissionSourceType = calculator::cell::TotalAggregatedFissionSourceMock;
  using MomentsType = system::moments::SphericalHarmonicMock;

  // Supporting objects
  std::unique_ptr<TotalAggregatedFissionSourceType> fission_source_mock_ptr_;
  std::unique_ptr<MomentsType> current_moments_ptr_;
  system::System test_system;

  void SetUp() override;
};

void EigKEffUpdaterViaFissionSourceTest::SetUp() {
  fission_source_mock_ptr_ =
      std::make_unique<TotalAggregatedFissionSourceType>();
  test_system.current_moments = std::move(current_moments_ptr_);
}

TEST_F(EigKEffUpdaterViaFissionSourceTest, Constructor) {
  eigenvalue::k_eigenvalue::UpdaterViaFissionSource
      test_k_eff_updater(std::move(fission_source_mock_ptr_), 1.0, 5.0);

  EXPECT_EQ(test_k_eff_updater.k_effective(), std::nullopt);
  EXPECT_EQ(test_k_eff_updater.initial_k_effective(), 1.0);
  EXPECT_EQ(test_k_eff_updater.initial_fission_source(), 5.0);
  EXPECT_EQ(test_k_eff_updater.current_fission_source(), std::nullopt);

  auto fission_source_obs_ptr_ =
      dynamic_cast<TotalAggregatedFissionSourceType*>(
          test_k_eff_updater.fission_source_calculator());

  EXPECT_NE(fission_source_obs_ptr_, nullptr);
}

TEST_F(EigKEffUpdaterViaFissionSourceTest, BadConstructor) {

  std::array<double,2> bad_keff{0, -1.3};
  std::array<double,2> bad_fission_sources{0, -10.0};

  for (auto const keff : bad_keff) {
    for (auto const fission_source : bad_fission_sources) {
      EXPECT_ANY_THROW({eigenvalue::k_eigenvalue::UpdaterViaFissionSource
            test_k_eff_updater(std::move(fission_source_mock_ptr_), keff, fission_source);});
    }
  }
}

TEST_F(EigKEffUpdaterViaFissionSourceTest, BadFissionSources) {
  eigenvalue::k_eigenvalue::UpdaterViaFissionSource
      test_k_eff_updater(std::move(fission_source_mock_ptr_), 1.0, 5.0);
  auto fission_source_obs_ptr_ =
      dynamic_cast<TotalAggregatedFissionSourceType*>(
          test_k_eff_updater.fission_source_calculator());

  EXPECT_CALL(*fission_source_obs_ptr_,
              AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(-1));

  EXPECT_ANY_THROW(test_k_eff_updater.CalculateK_Effective(test_system));

  EXPECT_CALL(*fission_source_obs_ptr_,
              AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(0));

  EXPECT_ANY_THROW(test_k_eff_updater.CalculateK_Effective(test_system));
}

TEST_F(EigKEffUpdaterViaFissionSourceTest, CalculateKEff) {
  eigenvalue::k_eigenvalue::UpdaterViaFissionSource
      test_k_eff_updater(std::move(fission_source_mock_ptr_), 1.5, 5.0);
  auto fission_source_obs_ptr_ =
      dynamic_cast<TotalAggregatedFissionSourceType*>(
          test_k_eff_updater.fission_source_calculator());

  const std::array<double, 3> fission_sources{5.3, 2.6, 1.3};

  EXPECT_CALL(*fission_source_obs_ptr_,
      AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(fission_sources[0]));

  double expected_result = 1.59;

  double k_effective = test_k_eff_updater.CalculateK_Effective(test_system);
  test_system.k_effective = k_effective;

  EXPECT_DOUBLE_EQ(k_effective, expected_result);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.current_fission_source().value_or(0),
      fission_sources[0]);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.k_effective().value_or(0),
      expected_result);

  EXPECT_CALL(*fission_source_obs_ptr_,
              AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(fission_sources[1]));

  k_effective = test_k_eff_updater.CalculateK_Effective(test_system);
  test_system.k_effective = k_effective;

  expected_result = 0.78;

  EXPECT_NEAR(k_effective, expected_result, 1e-12);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.current_fission_source().value_or(0),
      fission_sources[1]);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.k_effective().value_or(0),
      k_effective);

  EXPECT_CALL(*fission_source_obs_ptr_,
              AggregatedFissionSource(test_system.current_moments.get()))
      .WillOnce(Return(fission_sources[2]));

  k_effective = test_k_eff_updater.CalculateK_Effective(test_system);
  test_system.k_effective = k_effective;

  expected_result = 0.39;

  EXPECT_NEAR(k_effective, expected_result, 1e-12);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.current_fission_source().value_or(0),
      fission_sources[2]);
  EXPECT_DOUBLE_EQ(test_k_eff_updater.k_effective().value_or(0), k_effective);

}

} // namespace