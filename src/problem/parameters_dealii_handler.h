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
  ParametersDealiiHandler();
  ~ParametersDealiiHandler() = default;

  /*! \brief Parses a ParameterHandler object for problem parameters */
  void Parse(const dealii::ParameterHandler &handler);
  
  /*! \brief Set up a ParameterHandler object for problem parameters.
   * This setup includes default values.
   */
  void SetUp(dealii::ParameterHandler &handler);
  
  // Basic Parameters
  /*! Get problem spatial dimension */
  int  SpatialDimension() const override { return spatial_dimension_; }

 private:
  int spatial_dimension_;

  // Key-words for input file
  const std::string kSpatialDimension_ = "problem dimension";
};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
