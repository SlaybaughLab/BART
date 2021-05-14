#ifndef BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_HPP_
#define BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_HPP_

#include "test_helpers/gmock_wrapper.h"

#include "problem/parameters_i.hpp"

namespace bart {

namespace problem {

class ParametersMock : public ParametersI {
 public:
  MOCK_METHOD(int, NEnergyGroups, (), (const));
  MOCK_METHOD(EquationType, TransportModel, (), (const));

  MOCK_METHOD(bool, HaveReflectiveBC, (), (const));
  MOCK_METHOD(std::vector<int>, NCells, (), (const));
  MOCK_METHOD((std::map<Boundary, bool>), ReflectiveBoundary, (), (const));
  MOCK_METHOD(int, SpatialDimension, (), (const));
  MOCK_METHOD(std::vector<double>, SpatialMax, (), (const));
  MOCK_METHOD(int, UniformRefinements, (), (const));

  MOCK_METHOD(DiscretizationType, Discretization, (), (const));
  MOCK_METHOD(int, FEPolynomialDegree, (), (const));

  MOCK_METHOD(int, NumberOfMaterials, (), (const));
  MOCK_METHOD(std::string, MaterialMapFilename, (), (const));
  MOCK_METHOD((std::unordered_map<int, std::string>), MaterialFilenames, (), (const));

  MOCK_METHOD(bool, UseTwoGridAcceleration, (), (const, override));
  MOCK_METHOD(bool, DoNDA, (), (const));

  MOCK_METHOD(EigenSolverType, EigenSolver, (), (const));
  MOCK_METHOD(eigenvalue::k_eigenvalue::K_EffectiveUpdaterName, K_EffectiveUpdaterType, (), (const));
  MOCK_METHOD(InGroupSolverType, InGroupSolver, (), (const));
  MOCK_METHOD(LinearSolverType, LinearSolver, (), (const));

  MOCK_METHOD(AngularQuadType, AngularQuad, (), (const));
  MOCK_METHOD(int, AngularQuadOrder, (), (const));

  MOCK_METHOD(bool, DoDiscreteFourierTransformOfError, (), (const, override));
  MOCK_METHOD(std::string, OutputFilenameBase, (), (const));

  MOCK_METHOD(bool, OutputAggregatedSourceData, (), (const, override));
  MOCK_METHOD(bool, OutputScalarFluxAsVTU, (), (const, override));
  MOCK_METHOD(bool, OutputFissionSourceAsVTU, (), (const, override));
  MOCK_METHOD(bool, OutputScatteringSourceAsVTU, (), (const, override));
  MOCK_METHOD(bool, OutputInnerIterationsToFile, (), (const, override));

};

} // namespace problem 

} // namespace bart

#endif // BART_SRC_PROBLEM_TESTS_PARAMETERS_MOCK_HPP_
