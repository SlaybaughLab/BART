#include "../parameters_dealii_handler.h"

#include <stdexcept>
#include <vector>

#include "gtest/gtest.h"

#include "../parameter_types.h"

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_parameter_handler;
  bart::problem::ParametersDealiiHandler test_parameters;

  // Key-words for input file
  const std::string kNCells = "number of cells for x, y, z directions";
  const std::string kOutputFilenameBase = "output file name base";
  const std::string kSpatialDimension_ = "problem dimension";
  const std::string kSpatialMax_ = "x, y, z max values of boundary locations";
  const std::string kTransportModel_ = "transport model";
  
  const std::string kEigenSolver_ = "eigen solver name";
  const std::string kInGroupSolver_ = "in group solver name";
  const std::string kLinearSolver_ = "ho linear solver name";
  const std::string kMultiGroupSolver_ = "mg solver name";

  const std::string kAngularQuadrature_ = "angular quadrature name";
  const std::string kAngularQuadOrder_ = "angular quadrature order";
};

void ParametersDealiiHandlerTest::SetUp() {
  test_parameters.SetUp(test_parameter_handler);
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.SpatialDimension(), 2)
      << "Default spatial dimension";
  ASSERT_EQ(test_parameters.OutputFilenameBase(), "bart_output")
      << "Default spatial dimension";
  ASSERT_EQ(test_parameters.TransportModel(), bart::problem::EquationType::kNone)
      << "Default transport model";

  ASSERT_EQ(test_parameters.LinearSolver(),
            bart::problem::LinearSolverType::kConjugateGradient)
      << "Default linear solver";
  ASSERT_EQ(test_parameters.InGroupSolver(),
            bart::problem::InGroupSolverType::kSourceIteration)
      << "Default in-group solver";
  ASSERT_EQ(test_parameters.EigenSolver(),
            bart::problem::EigenSolverType::kPowerIteration)
      << "Default eigenvalue solver";
  ASSERT_EQ(test_parameters.MultiGroupSolver(),
            bart::problem::MultiGroupSolverType::kGaussSeidel)
      << "Default multi-group solver";
  

  ASSERT_EQ(test_parameters.AngularQuad(),
            bart::problem::AngularQuadType::kNone)
      << "Default angular quadrature";
  ASSERT_EQ(test_parameters.AngularQuadOrder(), 4)
      << "Default angular quadrature order";
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersParse) {

  std::vector<int>    n_cells{10, 5, 20};
  std::vector<double> spatial_max{10.0, 5.0, 8.0};
  std::string         output_filename_base{"test_output"};

  // Set testing Parameters
  test_parameter_handler.set(kNCells, "10, 5, 20");
  test_parameter_handler.set(kOutputFilenameBase, output_filename_base);
  test_parameter_handler.set(kSpatialDimension_, 3.0);
  test_parameter_handler.set(kSpatialMax_, "10.0, 5.0, 8.0");
  test_parameter_handler.set(kTransportModel_, "saaf");

  test_parameter_handler.set(kEigenSolver_, "none");
  test_parameter_handler.set(kInGroupSolver_, "none");
  test_parameter_handler.set(kLinearSolver_, "gmres");
  test_parameter_handler.set(kMultiGroupSolver_, "none");
  
  test_parameter_handler.set(kAngularQuadrature_, "lsgc");
  test_parameter_handler.set(kAngularQuadOrder_, "8");
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.NCells(), n_cells)
      << "Parsed number of cells";                             
  ASSERT_EQ(test_parameters.OutputFilenameBase(), output_filename_base)
      << "Parsed output filename base";
  ASSERT_EQ(test_parameters.SpatialDimension(), 3.0)
      << "Parsed spatial dimension";
  ASSERT_EQ(test_parameters.SpatialMax(), spatial_max)            
      << "Parsed spatial maximums";
  ASSERT_EQ(test_parameters.TransportModel(),
            bart::problem::EquationType::kSelfAdjointAngularFlux)
      << "Parsed transport model";

  ASSERT_EQ(test_parameters.EigenSolver(),
            bart::problem::EigenSolverType::kNone)
      << "Parsed eigenvalue solver";
    ASSERT_EQ(test_parameters.InGroupSolver(),
            bart::problem::InGroupSolverType::kNone)
      << "Parsed in-group solver";
    ASSERT_EQ(test_parameters.LinearSolver(),
            bart::problem::LinearSolverType::kGMRES)
      << "Parsed linear solver";
  ASSERT_EQ(test_parameters.MultiGroupSolver(),
            bart::problem::MultiGroupSolverType::kNone)
      << "Parsed multi-group solver";

  ASSERT_EQ(test_parameters.AngularQuad(),
            bart::problem::AngularQuadType::kLevelSymmetricGaussChebyshev)
      << "Parsed angular quadrature";
  ASSERT_EQ(test_parameters.AngularQuadOrder(), 8)
      << "Parsed angular quadrature order";
}
