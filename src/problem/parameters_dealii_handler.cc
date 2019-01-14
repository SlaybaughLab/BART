#include "parameters_dealii_handler.h"

#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {

  std::vector<double> double_vector;
  double read_value;
  
  spatial_dimension_ = handler.get_integer(kSpatialDimension_);
  // Parse string into vector<double>
  std::stringstream iss(handler.get(kSpatialMax_));
  std::cout << handler.get(kSpatialMax_);
  while(iss >> read_value) {
    double_vector.push_back(read_value);
    if (iss.peek() == ',' || iss.peek() == ' ')
      iss.ignore();
  }
  spatial_max = double_vector;
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  // Basic parameters
  handler.declare_entry(kSpatialDimension_, "2", Pattern::Integer(1, 3), "");
  handler.declare_entry (kSpatialMax_, "", Pattern::List (Pattern::Double ()),
        "xmax, ymax, zmax of the boundaries, mins are zero");
}

} // namespace problem

} // namespace bart
