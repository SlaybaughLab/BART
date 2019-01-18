#include "../parameters_dealii_handler.h"

#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

#include "../parameter_types.h"

class ParametersDealiiHandlerTest : public ::testing::Test {
 protected:
  void SetUp() override;
  dealii::ParameterHandler test_parameter_handler;
  bart::problem::ParametersDealiiHandler test_parameters;
  bart::problem::ParametersDealiiHandler::KeyWords key_words;
};

void ParametersDealiiHandlerTest::SetUp() {
  test_parameters.SetUp(test_parameter_handler);
}

TEST_F(ParametersDealiiHandlerTest, BasicParametersDefault) {
  std::map<bart::problem::Boundary, bool> test_reflective_map {
    {bart::problem::Boundary::kXMin, false},
    {bart::problem::Boundary::kXMax, false},
    {bart::problem::Boundary::kYMin, false},
    {bart::problem::Boundary::kYMax, false},
    {bart::problem::Boundary::kZMin, false},
    {bart::problem::Boundary::kZMax, false},
        };
  
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.Discretization(),
            bart::problem::DiscretizationType::kContinuousFEM)
      << "Default discretization";
  ASSERT_EQ(test_parameters.FEPolynomialDegree(), 1)
      << "Default finite element polynomial degree";
  ASSERT_EQ(test_parameters.IsEigenvalueProblem(), false)
      << "Default eigenvalue problem";
  ASSERT_EQ(test_parameters.HaveReflectiveBC(), false)
      << "Default have reflective boundaries";
  ASSERT_EQ(test_parameters.FirstThermalGroup(), 0)
      << "Default first thermal group";
  ASSERT_EQ(test_parameters.NEnergyGroups(), 1)
      << "Default number of energy groups";
  ASSERT_EQ(test_parameters.SpatialDimension(), 2)
      << "Default spatial dimension";
  ASSERT_EQ(test_parameters.ReflectiveBoundary(), test_reflective_map)
      << "Default first reflective boundaries";
  ASSERT_EQ(test_parameters.OutputFilenameBase(), "bart_output")
      << "Default spatial dimension";
  ASSERT_EQ(test_parameters.TransportModel(), bart::problem::EquationType::kNone)
      << "Default transport model";
}

TEST_F(ParametersDealiiHandlerTest, MeshParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.IsMeshGenerated(), true)
      << "Default mesh generation";
  ASSERT_EQ(test_parameters.MeshFilename(), "")
      << "Default mesh filename";
  ASSERT_EQ(test_parameters.UniformRefinements(), 0)
      << "Default number of uniform refinements";
  ASSERT_EQ(test_parameters.FuelRodRadius(), 0.5)
      << "Default fuel rod radius";
  ASSERT_EQ(test_parameters.FuelRodTriangulation(),
            bart::problem::FuelRodTriangulationType::kNone)
      << "Default fuel rod triangulation type";
}

TEST_F(ParametersDealiiHandlerTest, MaterialParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  std::unordered_map<int, std::string> empty_map;
  
  ASSERT_EQ(test_parameters.MaterialFilenames(), empty_map)
      << "Default material filenames";
  ASSERT_EQ(test_parameters.MaterialMapFilename(), "")
      << "Default mesh filename";
  ASSERT_EQ(test_parameters.NumberOfMaterials(), 1)
      << "Default number of materials";
}

TEST_F(ParametersDealiiHandlerTest, AccelerationParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.Preconditioner(),
            bart::problem::PreconditionerType::kAMG)
        << "Default preconditioner";
  ASSERT_EQ(test_parameters.BlockSSORFactor(), 1.0)
      << "Default BSSOR Factor"; 
  ASSERT_EQ(test_parameters.DoNDA(), false)
      << "Default NDA usage";
  ASSERT_EQ(test_parameters.NDADiscretization(),
            bart::problem::DiscretizationType::kContinuousFEM)
      << "Default NDA Discretization";
  ASSERT_EQ(test_parameters.NDALinearSolver(),
            bart::problem::LinearSolverType::kNone)
      << "Default NDA linear solver";
  ASSERT_EQ(test_parameters.NDAPreconditioner(),
            bart::problem::PreconditionerType::kJacobi)
        << "Default NDA preconditioner";
  ASSERT_EQ(test_parameters.NDABlockSSORFactor(), 1.0)
      << "Default NDA BSSOR Factor"; 
}

TEST_F(ParametersDealiiHandlerTest, SolverParametersDefault) {
  test_parameters.Parse(test_parameter_handler);

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

}

TEST_F(ParametersDealiiHandlerTest, AngularQuadParametersDefault) {
  test_parameters.Parse(test_parameter_handler);
  
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
  std::map<bart::problem::Boundary, bool> parsed_reflective_map {
    {bart::problem::Boundary::kXMin, true},
    {bart::problem::Boundary::kXMax, false},
    {bart::problem::Boundary::kYMin, false},
    {bart::problem::Boundary::kYMax, true},
    {bart::problem::Boundary::kZMin, false},
    {bart::problem::Boundary::kZMax, false},
        };

  // Set testing Parameters
  test_parameter_handler.set(key_words.kDiscretization_, "dfem");
  test_parameter_handler.set(key_words.kEigenvalueProblem_, "true");
  test_parameter_handler.set(key_words.kHaveReflectiveBC_, "true");
  test_parameter_handler.set(key_words.kFEPolynomialDegree_, "2");
  test_parameter_handler.set(key_words.kFirstThermalGroup_, "2");
  test_parameter_handler.set(key_words.kNCells_, "10, 5, 20");
  test_parameter_handler.set(key_words.kNEnergyGroups_, "10");
  test_parameter_handler.set(key_words.kOutputFilenameBase_, output_filename_base);
  test_parameter_handler.set(key_words.kReflectiveBoundary_, "xmin, ymax");
  test_parameter_handler.set(key_words.kSpatialDimension_, 3.0);
  test_parameter_handler.set(key_words.kSpatialMax_, "10.0, 5.0, 8.0");
  test_parameter_handler.set(key_words.kTransportModel_, "saaf");
  
  test_parameters.Parse(test_parameter_handler);


  ASSERT_EQ(test_parameters.Discretization(),
            bart::problem::DiscretizationType::kDiscontinuousFEM)
      << "Parsed discretization";
  ASSERT_EQ(test_parameters.IsEigenvalueProblem(), true)
      << "Parsed eigenvalue problem";
  ASSERT_EQ(test_parameters.FEPolynomialDegree(), 2)
        << "Parsed finite element polynomial degree";
  ASSERT_EQ(test_parameters.HaveReflectiveBC(), true)
      << "Parsed have reflective boundaries";
  ASSERT_EQ(test_parameters.FirstThermalGroup(), 2)
      << "Parsed first thermal group";
  ASSERT_EQ(test_parameters.NCells(), n_cells)
      << "Parsed number of cells";
  ASSERT_EQ(test_parameters.NEnergyGroups(), 10)
      << "Parsed number of enery groups";
  ASSERT_EQ(test_parameters.OutputFilenameBase(), output_filename_base)
      << "Parsed output filename base";
  ASSERT_EQ(test_parameters.ReflectiveBoundary(), parsed_reflective_map)
      << "Default first thermal group";
  ASSERT_EQ(test_parameters.SpatialDimension(), 3.0)
      << "Parsed spatial dimension";
  ASSERT_EQ(test_parameters.SpatialMax(), spatial_max)            
      << "Parsed spatial maximums";
  ASSERT_EQ(test_parameters.TransportModel(),
            bart::problem::EquationType::kSelfAdjointAngularFlux)
      << "Parsed transport model";
}

TEST_F(ParametersDealiiHandlerTest, MeshParametersParsed) {

  test_parameter_handler.set(key_words.kMeshGenerated_, "false");
  test_parameter_handler.set(key_words.kUniformRefinements_, "1");
  test_parameter_handler.set(key_words.kMeshFilename_, "test_mesh.msh");
  test_parameter_handler.set(key_words.kFuelRodRadius_, "1.0");
  test_parameter_handler.set(key_words.kFuelRodTriangulation_, "simple");

  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.IsMeshGenerated(), false)
      << "Parsed mesh generation";
  ASSERT_EQ(test_parameters.MeshFilename(), "test_mesh.msh")
      << "Parsed mesh file name";
  ASSERT_EQ(test_parameters.UniformRefinements(), 1)
      << "Parsed number of uniform refinements";
  ASSERT_EQ(test_parameters.FuelRodRadius(), 1.0)
      << "Default fuel rod radius";
  ASSERT_EQ(test_parameters.FuelRodTriangulation(),
            bart::problem::FuelRodTriangulationType::kSimple)
      << "Parsed fuel rod triangulation type";
}

TEST_F(ParametersDealiiHandlerTest, MaterialParametersParsed) {

  std::unordered_map<int, std::string> material_map{
    {1, "file_1"},
    {2, "file_2"}
  };
  
  test_parameter_handler.set(key_words.kNumberOfMaterials_, "5");
  
  test_parameter_handler.enter_subsection(key_words.kMaterialSubsection_);
  test_parameter_handler.set(key_words.kMaterialMapFilename_, "mid.txt");
  test_parameter_handler.set(key_words.kMaterialFilenames_,
                             "1: file_1, 2: file_2");
  test_parameter_handler.leave_subsection();
  
  test_parameters.Parse(test_parameter_handler);

  ASSERT_EQ(test_parameters.MaterialFilenames(), material_map)
      << "Parsed material filenames";
  ASSERT_EQ(test_parameters.MaterialMapFilename(), "mid.txt")
      << "Parsed material map file name";
  ASSERT_EQ(test_parameters.NumberOfMaterials(), 5)
      << "Parsed number of materials";                             
}

TEST_F(ParametersDealiiHandlerTest, AccelerationParametersParsed) {
  test_parameter_handler.set(key_words.kPreconditioner_, "bjacobi");
  test_parameter_handler.set(key_words.kBSSOR_Factor_, "1.5");
  test_parameter_handler.set(key_words.kDoNDA_, "true");
  test_parameter_handler.set(key_words.kNDA_Discretization_, "dfem");
  test_parameter_handler.set(key_words.kNDALinearSolver_, "gmres");
  test_parameter_handler.set(key_words.kNDAPreconditioner_, "amg");
  test_parameter_handler.set(key_words.kNDA_BSSOR_Factor_, "2.0");
  
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.Preconditioner(),
            bart::problem::PreconditionerType::kBlockJacobi)
        << "Parsed preconditioner";
  ASSERT_EQ(test_parameters.BlockSSORFactor(), 1.5)
      << "Parsed BSSOR Factor"; 
  ASSERT_EQ(test_parameters.DoNDA(), true)
      << "Parsed NDA usage";
  ASSERT_EQ(test_parameters.NDADiscretization(),
            bart::problem::DiscretizationType::kDiscontinuousFEM)
      << "Parsed NDA Discretization";
  ASSERT_EQ(test_parameters.NDALinearSolver(),
            bart::problem::LinearSolverType::kGMRES)
      << "Parsed NDA linear solver";
  ASSERT_EQ(test_parameters.NDAPreconditioner(),
            bart::problem::PreconditionerType::kAMG)
        << "Parsed NDA preconditioner";
  ASSERT_EQ(test_parameters.NDABlockSSORFactor(), 2.0)
      << "Parsed NDA BSSOR Factor"; 
}

TEST_F(ParametersDealiiHandlerTest, SolverParametersParsed) {

  test_parameter_handler.set(key_words.kEigenSolver_, "none");
  test_parameter_handler.set(key_words.kInGroupSolver_, "none");
  test_parameter_handler.set(key_words.kLinearSolver_, "gmres");
  test_parameter_handler.set(key_words.kMultiGroupSolver_, "none");
  
  test_parameters.Parse(test_parameter_handler);
  
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

}

TEST_F(ParametersDealiiHandlerTest, AngularQuadParametersParsed) {

  test_parameter_handler.set(key_words.kAngularQuad_, "lsgc");
  test_parameter_handler.set(key_words.kAngularQuadOrder_, "8");
  
  test_parameters.Parse(test_parameter_handler);
  
  ASSERT_EQ(test_parameters.AngularQuad(),
            bart::problem::AngularQuadType::kLevelSymmetricGaussChebyshev)
      << "Parsed angular quadrature";
  ASSERT_EQ(test_parameters.AngularQuadOrder(), 8)
      << "Parsed angular quadrature order";
}
