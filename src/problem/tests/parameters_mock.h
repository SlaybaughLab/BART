#ifndef BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_
#define BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_

#include "test_helpers/gmock_wrapper.h"

#include "problem/parameters_i.h"

namespace bart {

namespace problem {

class ParametersMock : public ParametersI {
 public:
  MOCK_METHOD(bool, DoDiscreteFourierTransformOfError, (), (const, override));

  MOCK_METHOD(DiscretizationType, Discretization, (), (const));

  MOCK_METHOD(bool, IsEigenvalueProblem, (), (const));

  MOCK_METHOD(int, FEPolynomialDegree, (), (const));

  MOCK_METHOD(int, FirstThermalGroup, (), (const));

  MOCK_METHOD(bool, HaveReflectiveBC, (), (const));

  MOCK_METHOD(std::vector<int>, NCells, (), (const));

  MOCK_METHOD(int, NEnergyGroups, (), (const));

  MOCK_METHOD(std::string, OutputFilenameBase, (), (const));

  MOCK_METHOD((std::map<Boundary, bool>), ReflectiveBoundary, (), (const));

  MOCK_METHOD(int, SpatialDimension, (), (const));

  MOCK_METHOD(std::vector<double>, SpatialMax, (), (const));

  MOCK_METHOD(EquationType, TransportModel, (), (const));
                                                                      
  MOCK_METHOD(bool, IsMeshGenerated, (), (const));

  MOCK_METHOD(std::string, MeshFilename, (), (const));

  MOCK_METHOD(int, UniformRefinements, (), (const));

  MOCK_METHOD(double, FuelPinRadius, (), (const));

  MOCK_METHOD(FuelPinTriangulationType, FuelPinTriangulation, (), (const));

  MOCK_METHOD(bool, IsMeshPinResolved, (), (const));

  MOCK_METHOD(int, NumberOfMaterials, (), (const));

  MOCK_METHOD(std::string, MaterialMapFilename, (), (const));

  MOCK_METHOD((std::unordered_map<int, std::string>), MaterialFilenames, (), (const));
  
  MOCK_METHOD(std::string, FuelPinMaterialMapFilename, (), (const));

  MOCK_METHOD(PreconditionerType, Preconditioner, (), (const));

  MOCK_METHOD(double, BlockSSORFactor, (), (const));

  MOCK_METHOD(bool, DoNDA, (), (const));

  MOCK_METHOD(DiscretizationType, NDADiscretization, (), (const));

  MOCK_METHOD(LinearSolverType, NDALinearSolver, (), (const));

  MOCK_METHOD(PreconditionerType, NDAPreconditioner, (), (const));

  MOCK_METHOD(double, NDABlockSSORFactor, (), (const));

  MOCK_METHOD(EigenSolverType, EigenSolver, (), (const));

  MOCK_METHOD(eigenvalue::k_eigenvalue::K_EffectiveUpdaterName, K_EffectiveUpdaterType, (), (const));

  MOCK_METHOD(InGroupSolverType, InGroupSolver, (), (const));

  MOCK_METHOD(LinearSolverType, LinearSolver, (), (const));

  MOCK_METHOD(MultiGroupSolverType, MultiGroupSolver, (), (const));

  MOCK_METHOD(AngularQuadType, AngularQuad, (), (const));

  MOCK_METHOD(int, AngularQuadOrder, (), (const));
  
};

} // namespace problem 

} // namespace bart

#endif // BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_
