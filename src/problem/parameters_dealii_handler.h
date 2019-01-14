#ifndef BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
#define BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_

#include <string>

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

class ParametersDealiiHandler : ParametersI {
 public:
  ParametersDealiiHandler(const dealii::ParameterHandler &handler);
  ~ParametersDealiiHandler() = default;

  // Basic Parameters
  int SpatialDimension() const override { return spatial_dimension_;} ;

 private:
  int spatial_dimension_;

  // Key-words for input file
  const std::string kw_dimension_ = "problem dimension";
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
