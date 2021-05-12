#include "problem/parameters_dealii_handler.hpp"

#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "problem/parameter_types.hpp"

namespace  {

using namespace bart::problem;

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  auto SetUp() -> void override;

  using K_EffectiveUpdaterName = bart::eigenvalue::k_eigenvalue::K_EffectiveUpdaterName;

  dealii::ParameterHandler test_parameter_handler{};
  ParametersDealiiHandler test_parameters{};
  ParametersDealiiHandler::KeyWords key_words{};
};

auto ParametersDealiiHandlerTest::SetUp() -> void {
  test_parameters.SetUp(test_parameter_handler);
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersDefault) {
  std::map<Boundary, bool> test_reflective_map {{Boundary::kXMin, false},{Boundary::kXMax, false},
                                                {Boundary::kYMin, false}, {Boundary::kYMax, false},
                                                {Boundary::kZMin, false}, {Boundary::kZMax, false}};
  
  test_parameters.Parse(test_parameter_handler);

  ASSERT_FALSE(test_parameters.DoDiscreteFourierTransformOfError()) << "Default do discrete fourier transform.";
  ASSERT_EQ(test_parameters.Discretization(), DiscretizationType::kContinuousFEM) << "Default discretization";
  ASSERT_EQ(test_parameters.FEPolynomialDegree(), 1) << "Default finite element polynomial degree";
  ASSERT_EQ(test_parameters.HaveReflectiveBC(), false) << "Default have reflective boundaries";
  ASSERT_EQ(test_parameters.NEnergyGroups(), 1) << "Default number of energy groups";
  ASSERT_EQ(test_parameters.SpatialDimension(), 2) << "Default spatial dimension";
  ASSERT_EQ(test_parameters.ReflectiveBoundary(), test_reflective_map) << "Default first reflective boundaries";
  ASSERT_EQ(test_parameters.OutputFilenameBase(), "bart_output") << "Default spatial dimension";
  ASSERT_EQ(test_parameters.TransportModel(), EquationType::kNone) << "Default transport model";
  EXPECT_FALSE(test_parameters.OutputAggregatedSourceData());
  EXPECT_FALSE(test_parameters.OutputScalarFluxAsVTU());
  EXPECT_FALSE(test_parameters.OutputFissionSourceAsVTU());
  EXPECT_FALSE(test_parameters.OutputScatteringSourceAsVTU());
  EXPECT_FALSE(test_parameters.OutputInnerIterationsToFile());
}

TEST_F(ParametersDealiiHandlerTest, MeshParametersDefault) {
  test_parameters.Parse(test_parameter_handler);
  ASSERT_EQ(test_parameters.UniformRefinements(), 0) << "Default number of uniform refinements";
}

TEST_F(ParametersDealiiHandlerTest, MaterialParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  std::unordered_map<int, std::string> empty_map;

  ASSERT_EQ(test_parameters.MaterialFilenames(), empty_map) << "Default material filenames";
  ASSERT_EQ(test_parameters.MaterialMapFilename(), "") << "Default mesh filename";
  ASSERT_EQ(test_parameters.NumberOfMaterials(), 1) << "Default number of materials";
}

TEST_F(ParametersDealiiHandlerTest, AccelerationParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.DoNDA(), false) << "Default NDA usage";
  ASSERT_EQ(test_parameters.UseTwoGridAcceleration(), false) << "Default two-grid usage";
}

TEST_F(ParametersDealiiHandlerTest, SolverParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.LinearSolver(), LinearSolverType::kGMRES) << "Default linear solver";
  ASSERT_EQ(test_parameters.InGroupSolver(), InGroupSolverType::kSourceIteration) << "Default in-group solver";
  ASSERT_EQ(test_parameters.EigenSolver(), EigenSolverType::kPowerIteration) << "Default eigenvalue solver";
  ASSERT_EQ(test_parameters.K_EffectiveUpdaterType(), K_EffectiveUpdaterName::kCalculatorViaFissionSource);
}

TEST_F(ParametersDealiiHandlerTest, AngularQuadParametersDefault) {
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.AngularQuad(), AngularQuadType::kNone) << "Default angular quadrature";
  ASSERT_EQ(test_parameters.AngularQuadOrder(), 4) << "Default angular quadrature order";
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersParse) {

  const std::vector<int>    n_cells{10, 5, 20};
  const std::vector<double> spatial_max{10.0, 5.0, 8.0};
  const std::string         output_filename_base{"test_output"};
  const std::map<Boundary, bool> parsed_reflective_map {{Boundary::kXMin, true}, {Boundary::kXMax, false},
                                                        {Boundary::kYMin, false}, {Boundary::kYMax, true},
                                                        {Boundary::kZMin, false}, {Boundary::kZMax, false}, };

  // Set testing Parameters
  test_parameter_handler.set(key_words.kDoDFTOfError_, "true");
  test_parameter_handler.set(key_words.kDiscretization_, "dfem");
  test_parameter_handler.set(key_words.kHaveReflectiveBC_, "true");
  test_parameter_handler.set(key_words.kFEPolynomialDegree_, "2");
  test_parameter_handler.set(key_words.kNCells_, "10, 5, 20");
  test_parameter_handler.set(key_words.kNEnergyGroups_, "10");
  test_parameter_handler.set(key_words.kOutputFilenameBase_, output_filename_base);
  test_parameter_handler.set(key_words.kReflectiveBoundary_, "xmin, ymax");
  test_parameter_handler.set(key_words.kSpatialDimension_, static_cast<double>(spatial_max.size()));
  test_parameter_handler.set(key_words.kSpatialMax_, "10.0, 5.0, 8.0");
  test_parameter_handler.set(key_words.kTransportModel_, "saaf");
  test_parameter_handler.set(key_words.kOutputAggregatedSourceData_, "true");
  test_parameter_handler.set(key_words.kOutputScalarFluxAsVTU_, "true");
  test_parameter_handler.set(key_words.kOutputFissionSourceAsVTU_, "true");
  test_parameter_handler.set(key_words.kOutputScatteringSourceAsVTU_, "true");
  test_parameter_handler.set(key_words.kOutputInnerIterationsToFile_, "true");
  
  test_parameters.Parse(test_parameter_handler);

  ASSERT_TRUE(test_parameters.DoDiscreteFourierTransformOfError());
  ASSERT_EQ(test_parameters.Discretization(), DiscretizationType::kDiscontinuousFEM) << "Parsed discretization";
  ASSERT_EQ(test_parameters.FEPolynomialDegree(), 2) << "Parsed finite element polynomial degree";
  ASSERT_EQ(test_parameters.HaveReflectiveBC(), true) << "Parsed have reflective boundaries";
  ASSERT_EQ(test_parameters.NCells(), n_cells) << "Parsed number of cells";
  ASSERT_EQ(test_parameters.NEnergyGroups(), 10) << "Parsed number of enery groups";
  ASSERT_EQ(test_parameters.OutputFilenameBase(), output_filename_base) << "Parsed output filename base";
  ASSERT_EQ(test_parameters.ReflectiveBoundary(), parsed_reflective_map) << "Default first thermal group";
  ASSERT_EQ(test_parameters.SpatialDimension(), 3.0) << "Parsed spatial dimension";
  ASSERT_EQ(test_parameters.SpatialMax(), spatial_max) << "Parsed spatial maximums";
  ASSERT_EQ(test_parameters.TransportModel(), EquationType::kSelfAdjointAngularFlux) << "Parsed transport model";
  EXPECT_TRUE(test_parameters.OutputAggregatedSourceData());
  EXPECT_TRUE(test_parameters.OutputScalarFluxAsVTU());
  EXPECT_TRUE(test_parameters.OutputFissionSourceAsVTU());
  EXPECT_TRUE(test_parameters.OutputScatteringSourceAsVTU());
  EXPECT_TRUE(test_parameters.OutputInnerIterationsToFile());
}

TEST_F(ParametersDealiiHandlerTest, MeshParametersParsed) {
  test_parameter_handler.set(key_words.kUniformRefinements_, "1");
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.UniformRefinements(), 1) << "Parsed number of uniform refinements";
}

TEST_F(ParametersDealiiHandlerTest, MaterialParametersParsed) {
  std::unordered_map<int, std::string> material_map{{1, "file_1"}, {2, "file_2"} };
  test_parameter_handler.set(key_words.kNumberOfMaterials_, "5");
  
  test_parameter_handler.enter_subsection(key_words.kMaterialSubsection_);
  test_parameter_handler.set(key_words.kMaterialMapFilename_, "mid.txt");
  test_parameter_handler.set(key_words.kMaterialFilenames_, "1: file_1, 2: file_2");
  test_parameter_handler.leave_subsection();
  
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.MaterialFilenames(), material_map) << "Parsed material filenames";
  ASSERT_EQ(test_parameters.MaterialMapFilename(), "mid.txt") << "Parsed material map file name";
  ASSERT_EQ(test_parameters.NumberOfMaterials(), 5) << "Parsed number of materials";
}

TEST_F(ParametersDealiiHandlerTest, AccelerationParametersParsed) {
  test_parameter_handler.set(key_words.kDoNDA_, "true");
  test_parameter_handler.set(key_words.kUseTwoGridAcceleration_, "true");
  test_parameters.Parse(test_parameter_handler);
  

  ASSERT_EQ(test_parameters.DoNDA(), true) << "Parsed NDA usage";
  ASSERT_EQ(test_parameters.UseTwoGridAcceleration(), true) << "Parsed two-grid usage";
}

TEST_F(ParametersDealiiHandlerTest, SolverParametersParsed) {

  test_parameter_handler.set(key_words.kEigenSolver_, "none");
  test_parameter_handler.set(key_words.kInGroupSolver_, "none");
  test_parameter_handler.set(key_words.kLinearSolver_, "gmres");
  test_parameter_handler.set(key_words.kK_EffectiveUpdaterType_, "rayleigh quotient");
  
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.EigenSolver(), EigenSolverType::kNone) << "Parsed eigenvalue solver";
  ASSERT_EQ(test_parameters.K_EffectiveUpdaterType(), K_EffectiveUpdaterName::kCalculatorViaRayleighQuotient);
  ASSERT_EQ(test_parameters.InGroupSolver(), InGroupSolverType::kNone) << "Parsed in-group solver";
  ASSERT_EQ(test_parameters.LinearSolver(), LinearSolverType::kGMRES) << "Parsed linear solver";
}

TEST_F(ParametersDealiiHandlerTest, AngularQuadParametersParsed) {
  test_parameter_handler.set(key_words.kAngularQuad_, "level_symmetric_gaussian");
  test_parameter_handler.set(key_words.kAngularQuadOrder_, "8");
  
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.AngularQuad(), AngularQuadType::kLevelSymmetricGaussian) << "Parsed angular quadrature";
  ASSERT_EQ(test_parameters.AngularQuadOrder(), 8) << "Parsed angular quadrature order";
}

} // namespace