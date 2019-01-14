#include "parameters_dealii_handler.h"

#include <deal.II/base/parameter_handler.h>

#include <unordered_map>
#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {

  std::map<std::string, EquationType> equation_type_map {
    {"ep",   EquationType::kEvenParity},
    {"saaf", EquationType::kSelfAdjointAngularFlux},
    {"none", EquationType::kNone}
  };
  
  n_cells_ = ParseDealiiIntList(handler.get(kNCells_));
  output_filename_base_ = handler.get(kOutputFilenameBase_);  
  spatial_dimension_ = handler.get_integer(kSpatialDimension_);
  spatial_max = ParseDealiiList(handler.get(kSpatialMax_));
  transport_model_ = equation_type_map[handler.get(kTransportModel_)];
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  /* Set parameters. A try and catch wrapper is required for some. For these,
   * we are using the build-in validity checking for the parameter, but the
   * default values are not valid themselves. Mostly this means that there is
   * no default value.
   */
  try {
    handler.declare_entry (kNCells_ , "",
                           Pattern::List (Pattern::Integer (0), 1, 3),
                           "Geometry is hyper rectangle defined by how many cells exist per direction");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}
  handler.declare_entry(kOutputFilenameBase_, "bart_output",Pattern::Anything(),
                         "name base of the output file");
  handler.declare_entry(kSpatialDimension_, "2", Pattern::Integer(1, 3), "");
  try {
    handler.declare_entry (kSpatialMax_, "",
                           Pattern::List(Pattern::Double(), 1, 3),
                           "xmax, ymax, zmax of the boundaries, mins are zero");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}
  handler.declare_entry(kTransportModel_, "none",
                        Pattern::Selection("ep|saaf|none"),
                        "valid names such as ep");
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

std::vector<int> ParametersDealiiHandler::ParseDealiiIntList(
    std::string to_parse) {
  
  std::vector<double> double_vector(ParseDealiiList(to_parse));
  std::vector<int> return_vector(double_vector.cbegin(), double_vector.cend());
  return return_vector;
}


} // namespace problem

} // namespace bart
