#include "framework/builder/framework_validator.h"

#include "problem/tests/parameters_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/reporter/tests/basic_reporter_mock.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return,
::testing::UnorderedElementsAreArray;

using Part = framework::builder::FrameworkPart;

class FrameworkBuilderFrameworkValidatorTest : public ::testing::Test {
 public:
  framework::builder::FrameworkValidator test_validator;
  NiceMock<problem::ParametersMock> mock_parameters;

  void SetUp() override;
};

void FrameworkBuilderFrameworkValidatorTest::SetUp() {
  ON_CALL(mock_parameters, IsEigenvalueProblem())
      .WillByDefault(Return(true));
  ON_CALL(mock_parameters, TransportModel())
      .WillByDefault(Return(problem::EquationType::kDiffusion));
  ON_CALL(mock_parameters, HaveReflectiveBC())
      .WillByDefault(Return(false));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestNonEigenvalue) {
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, IsEigenvalueProblem())
      .WillOnce(Return(false));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::ScatteringSourceUpdate}));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestEigenvalue) {
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, IsEigenvalueProblem())
      .WillOnce(DoDefault());

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::FissionSourceUpdate,
                                         Part::ScatteringSourceUpdate}));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestSAAFWithoutReflective) {
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, TransportModel())
      .Times(::testing::AnyNumber())
      .WillRepeatedly(Return(problem::EquationType::kSelfAdjointAngularFlux));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::FissionSourceUpdate,
                                         Part::ScatteringSourceUpdate}));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestSAAFWithReflective) {
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, TransportModel())
      .WillOnce(Return(problem::EquationType::kSelfAdjointAngularFlux));
  EXPECT_CALL(mock_parameters, HaveReflectiveBC())
      .WillOnce(Return(true));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::FissionSourceUpdate,
                                         Part::ScatteringSourceUpdate,
                                         Part::AngularSolutionStorage}));
}

} // namespace
