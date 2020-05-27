#include "framework/builder/framework_validator.h"

#include "problem/tests/parameters_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/reporter/tests/basic_reporter_mock.h"

namespace  {

using namespace bart;

using ::testing::UnorderedElementsAreArray, ::testing::Return;

using Part = framework::builder::FrameworkPart;

class FrameworkBuilderFrameworkValidatorTest : public ::testing::Test {
 public:
  framework::builder::FrameworkValidator test_validator;
  problem::ParametersMock mock_parameters;
};

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
      .WillOnce(Return(true));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::FissionSourceUpdate,
                                         Part::ScatteringSourceUpdate}));
}

} // namespace
