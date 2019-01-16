#ifndef BART_SRC_PROBLEM_PARAMETERS_I_H_
#define BART_SRC_PROBLEM_PARAMETERS_I_H_

#include <string>
#include <map>
#include <vector>

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

  // Basic Problem Parameters
  virtual int                      FirstThermalGroup()  const = 0;
  virtual bool                     HaveReflectiveBC()   const = 0;
  virtual std::vector<int>         NCells()             const = 0;
  virtual int                      NEnergyGroups()      const = 0;
  virtual int                      NumberOfMaterials()  const = 0;
  virtual std::string              OutputFilenameBase() const = 0;
  virtual std::map<Boundary, bool> ReflectiveBoundary() const = 0;
  virtual int                      SpatialDimension()   const = 0;
  virtual std::vector<double>      SpatialMax()         const = 0;
  virtual EquationType             TransportModel()     const = 0;
                                   
  // Acceleration parameters
  virtual PreconditionerType       Preconditioner()     const = 0;
  virtual bool                     DoNDA()              const = 0;
  virtual LinearSolverType         NDALinearSolver()    const = 0;
  virtual PreconditionerType       NDAPreconditioner()  const = 0;
                                     
  // Solver parameters             
  virtual EigenSolverType          EigenSolver()        const = 0;
  virtual InGroupSolverType        InGroupSolver()      const = 0;
  virtual LinearSolverType         LinearSolver()       const = 0;
  virtual MultiGroupSolverType     MultiGroupSolver()   const = 0;
                                   
  // Angular quadrature parameters
  virtual AngularQuadType          AngularQuad()        const = 0;
  virtual int                      AngularQuadOrder()   const = 0;
  
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_I_H_
