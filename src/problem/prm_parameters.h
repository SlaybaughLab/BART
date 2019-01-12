#ifndef BART_SRC_PROBLEM_PRM_PARAMETERS_H_
#define BART_SRC_PROBLEM_PRM_PARAMETERS_H_

#include <deal.II/base/parameter_handler.h>

#include "parameters_i.h"

namespace bart {

namespace problem {

/*!
 * \brief Problem parameters derived using a dealii ParameterHandler object.
 * 
 * The ParameterHandler can be given directly or a file that can be parsed into
 * a ParameterHandler object can be passed.
 */

class PrmParameters : ParametersI {
 public:
  PrmParameters(const dealii::ParameterHandler &handler);
  ~PrmParameters() = default;

  // Basic Parameters
  int Dimension() const override { return dimension_;} ;

 private:
  int dimension_;
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PRM_PARAMETERS_H_
