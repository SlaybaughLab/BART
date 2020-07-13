#include "framework/builder/framework_validator.h"

#include "problem/tests/parameters_mock.h"
#include "test_helpers/gmock_wrapper.h"
#include "utility/reporter/tests/basic_reporter_mock.h"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return,
::testing::UnorderedElementsAreArray, ::testing::IsEmpty,
::testing::HasSubstr, ::testing::ReturnRef, ::testing::A, ::testing::AllOf,
::testing::_;

using Part = framework::builder::FrameworkPart;

class FrameworkBuilderFrameworkValidatorTest : public ::testing::Test {
 public:
  framework::builder::FrameworkValidator test_validator;
  NiceMock<problem::ParametersMock> mock_parameters;
  NiceMock<utility::reporter::BasicReporterMock> mock_reporter;

  void SetUp() override;
};

void FrameworkBuilderFrameworkValidatorTest::SetUp() {
  using Boundary = problem::Boundary;
  ON_CALL(mock_parameters, IsEigenvalueProblem())
      .WillByDefault(Return(true));
  ON_CALL(mock_parameters, TransportModel())
      .WillByDefault(Return(problem::EquationType::kDiffusion));
  ON_CALL(mock_parameters, ReflectiveBoundary())
      .WillByDefault(Return(std::map<Boundary, bool>{
          {Boundary::kXMin, false},
          {Boundary::kXMax, false}}));
  ON_CALL(mock_reporter, Instream(A<utility::reporter::Color>()))
      .WillByDefault(ReturnRef(mock_reporter));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestNonEigenvalue) {
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, IsEigenvalueProblem())
      .WillOnce(Return(false));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::ScatteringSourceUpdate}));
  EXPECT_FALSE(test_validator.HasUnneededParts());
  EXPECT_EQ(test_validator.UnneededParts().size(), 0);
  EXPECT_TRUE(test_validator.Parts().empty());
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
  EXPECT_TRUE(test_validator.Parts().empty());
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
  EXPECT_TRUE(test_validator.Parts().empty());
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ParseTestSAAFWithReflective) {
  using Boundary = problem::Boundary;
  EXPECT_FALSE(test_validator.HasNeededParts());

  EXPECT_CALL(mock_parameters, TransportModel())
      .WillOnce(Return(problem::EquationType::kSelfAdjointAngularFlux));
  EXPECT_CALL(mock_parameters, ReflectiveBoundary())
      .WillOnce(Return(std::map<Boundary, bool>{
          {Boundary::kXMin, true},
          {Boundary::kXMax, false}}));

  test_validator.Parse(mock_parameters);

  EXPECT_TRUE(test_validator.HasNeededParts());
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::FissionSourceUpdate,
                                         Part::ScatteringSourceUpdate,
                                         Part::AngularSolutionStorage}));
  EXPECT_TRUE(test_validator.Parts().empty());
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, AddPart) {
  test_validator.Parse(mock_parameters);
  EXPECT_TRUE(test_validator.Parts().empty());
  test_validator
      .AddPart(Part::ScatteringSourceUpdate)
      .AddPart(Part::FissionSourceUpdate);
  EXPECT_THAT(test_validator.NeededParts(),
              UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                         Part::FissionSourceUpdate}));
  EXPECT_THAT(test_validator.Parts(),
              UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                         Part::FissionSourceUpdate}));
  EXPECT_FALSE(test_validator.HasUnneededParts());
  EXPECT_TRUE(test_validator.UnneededParts().empty());
  test_validator.AddPart(Part::AngularSolutionStorage);
  EXPECT_THAT(test_validator.Parts(),
              UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                         Part::FissionSourceUpdate,
                                         Part::AngularSolutionStorage}));
  EXPECT_TRUE(test_validator.HasUnneededParts());
  EXPECT_THAT(test_validator.UnneededParts(),
              UnorderedElementsAreArray({Part::AngularSolutionStorage}));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationMissing) {
  test_validator.Parse(mock_parameters);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate,
                                Part::FissionSourceUpdate};

  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Red))
      .Times(::testing::AtLeast(2));
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Reset))
      .Times(::testing::AtLeast(2));
  EXPECT_CALL(mock_reporter, Instream(A<const std::string&>()))
      .Times(::testing::AtLeast(2))
      .WillRepeatedly(ReturnRef(mock_reporter));
  EXPECT_CALL(mock_reporter, Report(_,utility::reporter::Color::Yellow));

  test_validator.ReportValidation(mock_reporter);
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationPresent) {
  test_validator.Parse(mock_parameters);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate,
                                Part::FissionSourceUpdate};
  EXPECT_EQ(test_validator.NeededParts(), expected_parts);
  test_validator
      .AddPart(Part::ScatteringSourceUpdate)
      .AddPart(Part::FissionSourceUpdate);

  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Green))
      .Times(::testing::AtLeast(2));
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Reset))
      .Times(::testing::AtLeast(2));
  EXPECT_CALL(mock_reporter, Instream(A<const std::string&>()))
      .Times(::testing::AtLeast(2))
      .WillRepeatedly(ReturnRef(mock_reporter));

  test_validator.ReportValidation(mock_reporter);
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationMixed) {
  test_validator.Parse(mock_parameters);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate,
                                Part::FissionSourceUpdate};

  test_validator
      .AddPart(Part::ScatteringSourceUpdate)
      .AddPart(Part::AngularSolutionStorage);
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Green));
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Yellow));
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Red));
  EXPECT_CALL(mock_reporter, Instream(A<const std::string&>()))
      .Times(::testing::AtLeast(3))
      .WillRepeatedly(ReturnRef(mock_reporter));
  EXPECT_CALL(mock_reporter, Instream(utility::reporter::Color::Reset))
      .Times(::testing::AtLeast(3));
  EXPECT_CALL(mock_reporter, Report(_,utility::reporter::Color::Yellow));

  test_validator.ReportValidation(mock_reporter);
}


} // namespace
