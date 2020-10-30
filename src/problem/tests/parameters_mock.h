#ifndef BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_
#define BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_

#include "test_helpers/gmock_wrapper.h"

#include "problem/parameters_i.h"

namespace bart {

namespace problem {

class ParametersMock : public ParametersI {
 public:
  MOCK_METHOD(bool, DoDiscreteFourierTransformOfError, (), (const, override));

  MOCK_CONST_METHOD0(Discretization, DiscretizationType());

  MOCK_CONST_METHOD0(IsEigenvalueProblem, bool());

  MOCK_CONST_METHOD0(FEPolynomialDegree, int());

  MOCK_CONST_METHOD0(FirstThermalGroup, int());

  MOCK_CONST_METHOD0(HaveReflectiveBC, bool());

  MOCK_CONST_METHOD0(NCells, std::vector<int>());

  MOCK_CONST_METHOD0(NEnergyGroups, int());

  MOCK_CONST_METHOD0(OutputFilenameBase, std::string());

  MOCK_CONST_METHOD0(ReflectiveBoundary, std::map<Boundary, bool>());

  MOCK_CONST_METHOD0(SpatialDimension, int());

  MOCK_CONST_METHOD0(SpatialMax, std::vector<double>());

  MOCK_CONST_METHOD0(TransportModel, EquationType());
                                                                      
  MOCK_CONST_METHOD0(IsMeshGenerated, bool());

  MOCK_CONST_METHOD0(MeshFilename, std::string());

  MOCK_CONST_METHOD0(UniformRefinements, int());

  MOCK_CONST_METHOD0(FuelPinRadius, double());

  MOCK_CONST_METHOD0(FuelPinTriangulation, FuelPinTriangulationType());

  MOCK_CONST_METHOD0(IsMeshPinResolved, bool());

  MOCK_CONST_METHOD0(NumberOfMaterials, int());

  MOCK_CONST_METHOD0(MaterialMapFilename, std::string());

  MOCK_CONST_METHOD0(MaterialFilenames, std::unordered_map<int, std::string>());
  
  MOCK_CONST_METHOD0(FuelPinMaterialMapFilename, std::string());

  MOCK_CONST_METHOD0(Preconditioner, PreconditionerType());

  MOCK_CONST_METHOD0(BlockSSORFactor, double());

  MOCK_CONST_METHOD0(DoNDA, bool());

  MOCK_CONST_METHOD0(NDADiscretization, DiscretizationType());

  MOCK_CONST_METHOD0(NDALinearSolver, LinearSolverType());

  MOCK_CONST_METHOD0(NDAPreconditioner, PreconditionerType());

  MOCK_CONST_METHOD0(NDABlockSSORFactor, double());

  MOCK_CONST_METHOD0(EigenSolver, EigenSolverType());

  MOCK_CONST_METHOD0(InGroupSolver, InGroupSolverType());

  MOCK_CONST_METHOD0(LinearSolver, LinearSolverType());

  MOCK_CONST_METHOD0(MultiGroupSolver, MultiGroupSolverType());

  MOCK_CONST_METHOD0(AngularQuad, AngularQuadType());

  MOCK_CONST_METHOD0(AngularQuadOrder, int());
  
};

} // namespace problem 

} // namespace bart

#endif // BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_H_
