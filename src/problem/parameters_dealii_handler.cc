#include "parameters_dealii_handler.h"

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {
  spatial_dimension_ = handler.get_integer(kSpatialDimension_);
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  // Basic parameters
  handler.declare_entry(kSpatialDimension_, "2", Pattern::Integer(1, 3), "");
}

} // namespace problem

} // namespace bart
