#ifndef BART_SRC_PROBLEM_PARAMETERS_I_H_
#define BART_SRC_PROBLEM_PARAMETERS_I_H_

#include <string>
#include <vector>

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
  virtual std::vector<int>    NCells()             const = 0;
  virtual std::string         OutputFilenameBase() const = 0;
  virtual int                 SpatialDimension()   const = 0;
  virtual std::vector<double> SpatialMax()         const = 0;
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_I_H_
