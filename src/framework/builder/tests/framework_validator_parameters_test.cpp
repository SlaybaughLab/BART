#include "framework/builder/framework_validator.hpp"

#include "framework/framework_parameters.hpp"
#include "instrumentation/tests/instrument_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "test_helpers/test_helper_functions.h"
#include "problem/parameter_types.h"
#include "utility/colors.hpp"

namespace  {

using namespace bart;

namespace builder = bart::framework::builder;
namespace helper = bart::test_helpers;

using ::testing::NiceMock, ::testing::Pair, ::testing::_;

class FrameworkBuilderFrameworkValidatorParametersTest : public ::testing::Test {
 public:
  using Part = builder::FrameworkPart;
  using FrameworkValidator = builder::FrameworkValidator;
  using FrameworkParameters = bart::framework::FrameworkParameters;
  using StatusInstrumentType = instrumentation::InstrumentMock<std::pair<std::string, utility::Color>>;

  FrameworkParameters framework_parameters_;
  FrameworkValidator test_validator_;

  std::shared_ptr<NiceMock<StatusInstrumentType>> mock_instrument;

  auto SetUp() -> void override;
};

auto FrameworkBuilderFrameworkValidatorParametersTest::SetUp() -> void {
  framework_parameters_ = {.neutron_energy_groups = 2};
  mock_instrument = std::make_shared<NiceMock<StatusInstrumentType>>();
  using ValidationStatusPort = builder::data_port::ValidatorStatusPort;
  instrumentation::GetPort<ValidationStatusPort>(test_validator_).AddInstrument(mock_instrument);
}

TEST_F(FrameworkBuilderFrameworkValidatorParametersTest, DefaultFramework) {
  test_validator_.Parse(framework_parameters_);
  EXPECT_TRUE(test_validator_.HasNeededParts());
  EXPECT_TRUE(test_validator_.NeededParts().contains(Part::ScatteringSourceUpdate));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::FissionSourceUpdate));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::AngularSolutionStorage));
}

TEST_F(FrameworkBuilderFrameworkValidatorParametersTest, DefaultEigensolve) {
  framework_parameters_.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  test_validator_.Parse(framework_parameters_);
  EXPECT_TRUE(test_validator_.HasNeededParts());
  EXPECT_TRUE(test_validator_.NeededParts().contains(Part::ScatteringSourceUpdate));
  EXPECT_TRUE(test_validator_.NeededParts().contains(Part::FissionSourceUpdate));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::AngularSolutionStorage));
}

TEST_F(FrameworkBuilderFrameworkValidatorParametersTest, NoneTypeEigensolverSendsWarning) {
  framework_parameters_.eigen_solver_type = problem::EigenSolverType::kNone;
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kYellow)));
  test_validator_.Parse(framework_parameters_);
  EXPECT_TRUE(test_validator_.HasNeededParts());
  EXPECT_TRUE(test_validator_.NeededParts().contains(Part::ScatteringSourceUpdate));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::FissionSourceUpdate));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::AngularSolutionStorage));
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

TEST_F(FrameworkBuilderFrameworkValidatorParametersTest, SAAFRequiresAngularSolutionStorage) {
  FrameworkParameters saaf_parameters{framework_parameters_};
  saaf_parameters.equation_type = problem::EquationType::kSelfAdjointAngularFlux;
  EXPECT_NO_THROW(test_validator_.Parse(saaf_parameters));
  EXPECT_FALSE(test_validator_.NeededParts().contains(Part::AngularSolutionStorage));

  saaf_parameters.reflective_boundaries = {problem::Boundary::kXMin};
  EXPECT_NO_THROW(test_validator_.Parse(saaf_parameters));
  EXPECT_TRUE(test_validator_.NeededParts().contains(Part::AngularSolutionStorage));
}



} // namespace


