#include "instrumentation/converter/system/group_scalar_flux_extractor.hpp"

#include <array>

#include "system/moments/tests/spherical_harmonic_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace instrumentation = bart::instrumentation;
namespace moments = bart::system::moments;
namespace test_helpers = bart::test_helpers;

using ::testing::DoDefault, ::testing::Return, ::testing::ReturnRef;

class InstrumentationConverterGroupScalarFluxExtractorTest : public ::testing::Test {
 public:
  using DealiiVector = dealii::Vector<double>;
  using TestConverter = instrumentation::converter::system::GroupScalarFluxExtractor;
  using SphericalHarmonic = moments::SphericalHarmonicMock;

  std::unique_ptr<moments::SphericalHarmonicI> test_spherical_harmonics_;
  SphericalHarmonic* spherical_harmonics_obs_ptr_{ nullptr };

  const int total_groups_{ test_helpers::RandomInt(5, 10) };
  const int max_l_{ test_helpers::RandomInt(0, 5) };
  std::vector<DealiiVector> group_scalar_fluxes_{};

  void SetUp() override;
};

void InstrumentationConverterGroupScalarFluxExtractorTest::SetUp() {
  const int solution_size { test_helpers::RandomInt(5, 10) };
  test_spherical_harmonics_ = std::make_unique<SphericalHarmonic>();
  spherical_harmonics_obs_ptr_ =
      dynamic_cast<SphericalHarmonic*>(test_spherical_harmonics_.get());

  ON_CALL(*spherical_harmonics_obs_ptr_,
          total_groups()).WillByDefault(Return(total_groups_));

  for (int i = 0; i < total_groups_; ++i) {
    DealiiVector group_scalar_flux(solution_size);
    for (DealiiVector::size_type j = 0; j < group_scalar_flux.size(); ++j) {
      group_scalar_flux[j] = test_helpers::RandomDouble(-100, 100);
    }
    group_scalar_fluxes_.emplace_back(group_scalar_flux);
  }

  for (std::size_t group = 0; group < group_scalar_fluxes_.size(); ++group) {
    ON_CALL(*spherical_harmonics_obs_ptr_,
            GetMoment(std::array<int,3>{static_cast<int>(group), 0, 0}))
        .WillByDefault(ReturnRef(group_scalar_fluxes_.at(group)));
  }
}

TEST_F(InstrumentationConverterGroupScalarFluxExtractorTest,
       ConstructorAndGroupToExtractGetter) {
  const int group_to_extract{ test_helpers::RandomInt(0, 20) };
  TestConverter test_converter(group_to_extract);
  EXPECT_EQ(test_converter.group_to_extract(), group_to_extract);
}

TEST_F(InstrumentationConverterGroupScalarFluxExtractorTest, Convert) {
  const int group_to_extract{ test_helpers::RandomInt(0, total_groups_) };
  TestConverter test_converter(group_to_extract);
  EXPECT_CALL(*spherical_harmonics_obs_ptr_,
              GetMoment(std::array{group_to_extract, 0, 0}))
      .WillOnce(DoDefault());
  auto extracted_value = test_converter.Convert(*test_spherical_harmonics_);
  ASSERT_NE(extracted_value.size(), 0);
  EXPECT_EQ(extracted_value, group_scalar_fluxes_.at(group_to_extract));
}



} // namespace
