#ifndef BART_SRC_PROBLEM_PARAMETERS_I_H_
#define BART_SRC_PROBLEM_PARAMETERS_I_H_

#include <string>
#include <map>
#include <unordered_map>
#include <vector>

#include "eigenvalue/k_effective/factory.hpp"
#include "parameter_types.h"

namespace bart {

namespace problem {

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

  // Fourier analysis
  virtual bool DoDiscreteFourierTransformOfError() const = 0;

  // Basic Problem Parameters
  /*! \brief Gets the type of spatial discretization of the transport equation */
  virtual DiscretizationType         Discretization()                 const = 0;
  /*! \brief Gets if solving an eigenvalue problem */
  virtual bool                       IsEigenvalueProblem()            const = 0;
  /*! \brief Gets polynomial degree of the finite element basis */
  virtual int                        FEPolynomialDegree()             const = 0;
  /*! \brief Gets the first thermal group in the group structure */
  virtual int                        FirstThermalGroup()              const = 0;
  /*! \brief Gets if the problem has reflective boundary conditions */
  virtual bool                       HaveReflectiveBC()               const = 0;
  /*! \brief Gets the number of cells in the spatial discretization in each
   * dimension. */
  virtual std::vector<int>           NCells()                         const = 0;
  /*! \brief Gets the total number of energy groups in the problem */
  virtual int                        NEnergyGroups()                  const = 0;
  /*! \brief Gets the base filename used for generating output files */
  virtual std::string                OutputFilenameBase()             const = 0;
  /*! \brief Gets a mapping of each boundary to a bool indicating if it is
   * a reflective boundary */
  virtual std::map<Boundary, bool>   ReflectiveBoundary()             const = 0;
  /*! \brief Gets the number of spatial dimensions in the problem */
  virtual int                        SpatialDimension()               const = 0;
  /*! \brief Gets the size of the spatial dimensions */
  virtual std::vector<double>        SpatialMax()                     const = 0;
  /*! \brief Gets the form of the transport equation to solve */
  virtual EquationType               TransportModel()                 const = 0;
                                                                      
  // Mesh parameters
  /*! \brief Gets if bart will be generating a mesh using deal.II */
  virtual bool                       IsMeshGenerated()                const = 0;
  /*! \brief Gets the filename of the mesh if not generated via deal.II */
  virtual std::string                MeshFilename()                   const = 0;
  /*! \brief Gets the number of uniform refinements */
  virtual int                        UniformRefinements()             const = 0;
  /*! \brief Gets the radius of the fuel pin if the problem has them */
  virtual double                     FuelPinRadius()                  const = 0;
  /*! \brief Gets the triangulation type of the fuel pin if present */
  virtual FuelPinTriangulationType   FuelPinTriangulation()           const = 0;
  /*! \brief Gets if the problem is pin resolved */
  virtual bool                       IsMeshPinResolved()              const = 0;
                                                                      
  // Material parameters
  /*! \brief Gets total number of materials in the problem */
  virtual int                        NumberOfMaterials()              const = 0;
  /*! \brief Gets the filename that shows the layout of materials in the problem */
  virtual std::string                MaterialMapFilename()            const = 0;
  /*! \brief Gets the mapping of material ids to the filename that contains
   * their material properties */
  virtual std::unordered_map<int, std::string>                       
                                     MaterialFilenames()              const = 0;
  /*! \brief Gets the filename that shows the layout of materials in the fuel
   * pin if present */
  virtual std::string                FuelPinMaterialMapFilename()        const = 0;
                                                            
  // Acceleration parameters
  /*! \brief Gets the type of preconditioner to use */
  virtual PreconditionerType         Preconditioner()                 const = 0;
  /*! \brief Gets the damping factor for block SSOR if used */
  virtual double                     BlockSSORFactor()                const = 0;
  /*! \brief Gets if NDA should be used */
  virtual bool                       DoNDA()                          const = 0;
  /*! \brief Gets the NDA discretization to use */
  virtual DiscretizationType         NDADiscretization()              const = 0;
  /*! \brief Gets which solver to use for NDA */
  virtual LinearSolverType           NDALinearSolver()                const = 0;
  /*! \brief Gets which preconditioner to use for NDA */
  virtual PreconditionerType         NDAPreconditioner()              const = 0;
  /*! \brief Gets damping factor for block SSOR if used for NDA */
  virtual double                     NDABlockSSORFactor()             const = 0;
                                                                      
  // Solver parameters
  /*! \brief Gets solver type for eigen iterations */
  virtual EigenSolverType            EigenSolver()                    const = 0;
  /*! \brief Get type of k-effective updater */
  virtual auto K_EffectiveUpdaterType() const -> eigenvalue::k_effective::K_EffectiveUpdaterName = 0;
  /*! \brief Gets solver type for in-group solves */
  virtual InGroupSolverType          InGroupSolver()                  const = 0;
  /*! \brief Gets solver type for linear solves */
  virtual LinearSolverType           LinearSolver()                   const = 0;
  /*! \brief Gets solver type for multi-group solves */
  virtual MultiGroupSolverType       MultiGroupSolver()               const = 0;
                                                                      
  // Angular quadrature parameters
  /*! \brief Gets type of angular quadrature to use */
  virtual AngularQuadType            AngularQuad()                    const = 0;
  /*! \brief Gets angular quadrature order */
  virtual int                        AngularQuadOrder()               const = 0;
  
};



} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_I_H_
