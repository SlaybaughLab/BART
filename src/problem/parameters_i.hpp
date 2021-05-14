#ifndef BART_SRC_PROBLEM_PARAMETERS_I_HPP_
#define BART_SRC_PROBLEM_PARAMETERS_I_HPP_

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

#include "eigenvalue/k_eigenvalue/factory.hpp"
#include "parameter_types.hpp"

namespace bart::problem {

/*!
 * \brief Interface for an object that contains problem parameters.
 *
 * This translates an input (file, etc) into values used by BART to determine
 * how to build the problem: what transport model to use, what angular
 * quadrature, etc.
 * 
 */

class ParametersI {
 public:
  
  virtual ~ParametersI() = default;

  // Equation parameters
  /*! \brief Gets the total number of energy groups in the problem */
  virtual auto NEnergyGroups() const -> int = 0;
  /*! \brief Gets the form of the transport equation to solve */
  virtual auto TransportModel() const -> EquationType = 0;
                                                                      
  // Domain parameters
  /*! \brief Gets if the problem has reflective boundary conditions */
  virtual auto HaveReflectiveBC() const -> bool = 0;
  /*! \brief Gets the number of cells in the spatial discretization in each
   * dimension. */
  virtual auto NCells() const -> std::vector<int> = 0;
  /*! \brief Gets a mapping of each boundary to a bool indicating if it is
   * a reflective boundary */
  virtual auto ReflectiveBoundary() const -> std::map<Boundary, bool> = 0;
  /*! \brief Gets the number of spatial dimensions in the problem */
  virtual auto SpatialDimension() const -> int = 0;
  /*! \brief Gets the size of the spatial dimensions */
  virtual auto SpatialMax() const -> std::vector<double> = 0;
  /*! \brief Gets the number of uniform refinements */
  virtual auto UniformRefinements() const -> int = 0;

  // FEM Parameters
  /*! \brief Gets the type of spatial discretization of the transport equation */
  virtual auto Discretization() const -> DiscretizationType = 0;
  /*! \brief Gets polynomial degree of the finite element basis */
  virtual auto FEPolynomialDegree() const -> int = 0;

  // Material parameters
  /*! \brief Gets total number of materials in the problem */
  virtual auto NumberOfMaterials() const -> int = 0;
  /*! \brief Gets the filename that shows the layout of materials in the problem */
  virtual auto MaterialMapFilename() const -> std::string = 0;
  /*! \brief Gets the mapping of material ids to the filename that contains
   * their material properties */
  virtual auto MaterialFilenames() const -> std::unordered_map<int, std::string> = 0;
                                                            
  // Acceleration parameters
  /*! \brief Use two-grid acceleration. */
  virtual auto UseTwoGridAcceleration() const -> bool = 0;
  /*! \brief Gets if NDA should be used */
  virtual auto DoNDA() const -> bool = 0;
                                                                      
  // Solver parameters
  /*! \brief Gets solver type for eigen iterations */
  virtual auto EigenSolver() const -> EigenSolverType = 0;
  /*! \brief Get type of k-effective updater */
  virtual auto K_EffectiveUpdaterType() const -> eigenvalue::k_eigenvalue::K_EffectiveUpdaterName = 0;
  /*! \brief Gets solver type for in-group solves */
  virtual auto InGroupSolver() const -> InGroupSolverType = 0;
  /*! \brief Gets solver type for linear solves */
  virtual auto LinearSolver() const -> LinearSolverType = 0;
                                                                      
  // Quadrature
  /*! \brief Gets type of angular quadrature to use */
  virtual auto AngularQuad() const -> AngularQuadType = 0;
  /*! \brief Gets angular quadrature order */
  virtual auto AngularQuadOrder() const -> int = 0;

  // Instrumentation and output
  /*! \brief Gets the base filename used for generating output files */
  virtual auto DoDiscreteFourierTransformOfError() const -> bool = 0;
  virtual auto OutputFilenameBase() const -> std::string = 0;
  virtual auto OutputAggregatedSourceData() const -> bool = 0;
  virtual auto OutputScalarFluxAsVTU() const -> bool = 0;
  virtual auto OutputFissionSourceAsVTU() const -> bool = 0;
  virtual auto OutputScatteringSourceAsVTU() const -> bool = 0;
  virtual auto OutputInnerIterationsToFile() const -> bool = 0;
};

} // namespace bart::problem

#endif // BART_SRC_PROBLEM_PARAMETERS_I_HPP_
