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
  PrmParameters(dealii::ParameterHandler handler);
  ~PrmParameters() = default;
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PRM_PARAMETERS_H_
