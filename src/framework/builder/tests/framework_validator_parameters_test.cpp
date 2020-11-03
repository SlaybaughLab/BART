#include "framework/builder/framework_validator.hpp"

#include "framework/framework_parameters.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"

namespace  {

namespace builder = bart::framework::builder;
namespace helper = bart::test_helpers;

class FrameworkBuilderFrameworkValidatorParametersTest : public ::testing::Test {
 public:
  using Part = builder::FrameworkPart;
  using FrameworkValidator = builder::FrameworkValidator;
  using FrameworkParameters = bart::framework::FrameworkParameters;

  FrameworkParameters framework_parameters_;
  FrameworkValidator test_validator_;
  auto SetUp() -> void override;
};

auto FrameworkBuilderFrameworkValidatorParametersTest::SetUp() -> void {
  framework_parameters_ = {.neutron_energy_groups = 2};
}

TEST_F(FrameworkBuilderFrameworkValidatorParametersTest, BadEnergyGroups) {
  FrameworkParameters bad_parameters{framework_parameters_};

  std::vector<double> bad_groups{helper::RandomVector(10, -10, 0)};
  bad_groups.emplace_back(0);

  for (const int bad_group : bad_groups) {
    bad_parameters.neutron_energy_groups = helper::RandomInt(-20, -1);
    EXPECT_ANY_THROW(test_validator_.Parse(bad_parameters))
              << "No throw for neutron_energy_groups = " + std::to_string(bad_group);
  }
}



} // namespace


