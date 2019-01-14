#include "parameters_dealii_handler.h"

#include <deal.II/base/parameter_handler.h>

#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {


  output_filename_base_ = handler.get(kOutputFilenameBase);  
  spatial_dimension_ = handler.get_integer(kSpatialDimension_);
  spatial_max = ParseDealiiList(handler.get(kSpatialMax_));
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  // Basic parameters
  handler.declare_entry (kOutputFilenameBase, "bart_output",Pattern::Anything(),
                         "name base of the output file");
  handler.declare_entry(kSpatialDimension_, "2", Pattern::Integer(1, 3), "");

  // Entries that have a default value that does not match the required
  // pattern defintion (or no default value)
  try {
    handler.declare_entry (kSpatialMax_, "",
                           Pattern::List(Pattern::Double(), 1, 3),
                           "xmax, ymax, zmax of the boundaries, mins are zero");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}
}

std::vector<double> ParametersDealiiHandler::ParseDealiiList(
    std::string to_parse) {

  std::vector<double> return_vector;
  double read_value;
  std::stringstream iss{to_parse};
  
  while(iss >> read_value) {
    return_vector.push_back(read_value);
    // Ignore next if comma or space
    if (iss.peek() == ',' || iss.peek() == ' ')
      iss.ignore();
  }

  return return_vector;
}

} // namespace problem

} // namespace bart
