#include "framework/builder/framework_validator.hpp"

#include "framework/framework_parameters.hpp"
#include "test_helpers/gmock_wrapper.h"
#include "instrumentation/tests/instrument_mock.h"
#include "utility/colors.hpp"

namespace  {

using namespace bart;

using ::testing::DoDefault, ::testing::NiceMock, ::testing::Return,
::testing::UnorderedElementsAreArray, ::testing::IsEmpty,
::testing::HasSubstr, ::testing::ReturnRef, ::testing::A, ::testing::AllOf,
::testing::_, ::testing::ElementsAre, ::testing::Pair;

using Part = framework::builder::FrameworkPart;
//TODO: If parameter handler support is pulled from Parse, rename these tests Base tests
class FrameworkBuilderFrameworkValidatorTest : public ::testing::Test {
 public:
  using StatusInstrumentType = instrumentation::InstrumentMock<std::pair<std::string, utility::Color>>;
  std::unique_ptr<framework::builder::FrameworkValidator> test_validator;
  framework::FrameworkParameters framework_parameters_;
  std::shared_ptr<NiceMock<StatusInstrumentType>> mock_instrument;

  void SetUp() override;
};

void FrameworkBuilderFrameworkValidatorTest::SetUp() {
  using Boundary = problem::Boundary;
  framework_parameters_.eigen_solver_type = problem::EigenSolverType::kPowerIteration;
  framework_parameters_.equation_type = problem::EquationType::kDiffusion;
  framework_parameters_.reflective_boundaries = {Boundary::kXMin, Boundary::kXMax};

  test_validator = std::make_unique<framework::builder::FrameworkValidator>();
  mock_instrument = std::make_shared<NiceMock<StatusInstrumentType>>();
  using ValidationStatusPort = framework::builder::data_port::ValidatorStatusPort;
  instrumentation::GetPort<ValidationStatusPort>(*test_validator).AddInstrument(mock_instrument);
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, AddPart) {
  test_validator->Parse(framework_parameters_);

  EXPECT_TRUE(test_validator->Parts().empty());
  test_validator->AddPart(Part::ScatteringSourceUpdate).AddPart(Part::FissionSourceUpdate);

  EXPECT_THAT(test_validator->NeededParts(), UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                                                        Part::FissionSourceUpdate}));
  EXPECT_THAT(test_validator->Parts(), UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                                                  Part::FissionSourceUpdate}));
  EXPECT_FALSE(test_validator->HasUnneededParts());
  EXPECT_TRUE(test_validator->UnneededParts().empty());
  test_validator->AddPart(Part::AngularSolutionStorage);
  EXPECT_THAT(test_validator->Parts(), UnorderedElementsAreArray({Part::ScatteringSourceUpdate,
                                                                  Part::FissionSourceUpdate,
                                                                  Part::AngularSolutionStorage}));
  EXPECT_TRUE(test_validator->HasUnneededParts());
  EXPECT_THAT(test_validator->UnneededParts(), UnorderedElementsAreArray({Part::AngularSolutionStorage}));
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationMissing) {
  test_validator->Parse(framework_parameters_);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate, Part::FissionSourceUpdate};
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kReset))).Times(::testing::AtLeast(1));
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kRed))).Times(::testing::AtLeast(2));
  EXPECT_CALL(*mock_instrument, Read(Pair(_,utility::Color::kYellow)));

  test_validator->ReportValidation();
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationPresent) {
  test_validator->Parse(framework_parameters_);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate, Part::FissionSourceUpdate};
  EXPECT_EQ(test_validator->NeededParts(), expected_parts);
  test_validator->AddPart(Part::ScatteringSourceUpdate).AddPart(Part::FissionSourceUpdate);
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kReset))).Times(::testing::AtLeast(1));
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kGreen))).Times(::testing::AtLeast(2));
  test_validator->ReportValidation();
}

TEST_F(FrameworkBuilderFrameworkValidatorTest, ReportValidationMixed) {
  test_validator->Parse(framework_parameters_);

  std::set<Part> expected_parts{Part::ScatteringSourceUpdate, Part::FissionSourceUpdate};

  test_validator->AddPart(Part::ScatteringSourceUpdate).AddPart(Part::AngularSolutionStorage);
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kReset))).Times(::testing::AtLeast(1));
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kGreen)));
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kYellow))).Times(2);
  EXPECT_CALL(*mock_instrument, Read(Pair(_, utility::Color::kRed)));

  test_validator->ReportValidation();
}


} // namespace
