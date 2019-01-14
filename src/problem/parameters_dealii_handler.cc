#include "parameters_dealii_handler.h"

#include <deal.II/base/parameter_handler.h>

#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {

  // Parse parameters
  // Basic Parameters
  n_cells_ = ParseDealiiIntList(handler.get(kNCells_));
  output_filename_base_ = handler.get(kOutputFilenameBase_);  
  spatial_dimension_ = handler.get_integer(kSpatialDimension_);
  spatial_max = ParseDealiiList(handler.get(kSpatialMax_));
  transport_model_ = kEquationTypeMap_.at(handler.get(kTransportModel_));

  // Solvers
  eigen_solver_ = kEigenSolverTypeMap_.at(handler.get(kEigenSolver_));
  linear_solver_ = kLinearSolverTypeMap_.at(handler.get(kLinearSolver_));
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  /* Set parameters. A try and catch wrapper is required for some. For these,
   * we are using the build-in validity checking for the parameter, but the
   * default values are not valid themselves. Mostly this means that there is
   * no default value.
   */
  SetUpBasicParameters(handler);
  SetUpSolverParameters(handler);
}

void ParametersDealiiHandler::SetUpBasicParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  try {
    handler.declare_entry(kNCells_ , "",
                          Pattern::List (Pattern::Integer (0), 1, 3),
                          "Geometry is hyper rectangle defined by how many cells exist per direction");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}
  
  handler.declare_entry(kOutputFilenameBase_, "bart_output",Pattern::Anything(),
                         "name base of the output file");
  
  handler.declare_entry(kSpatialDimension_, "2", Pattern::Integer(1, 3), "");
  
  try {
    handler.declare_entry(kSpatialMax_, "",
                          Pattern::List(Pattern::Double(), 1, 3),
                          "xmax, ymax, zmax of the boundaries, mins are zero");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}

  std::string equation_options{GetOptionString(kEquationTypeMap_)};
  handler.declare_entry(kTransportModel_, "none",
                        Pattern::Selection(equation_options),
                        "valid names such as ep");
}

void ParametersDealiiHandler::SetUpSolverParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  std::string linear_solver_options{GetOptionString(kLinearSolverTypeMap_)};
  handler.declare_entry(kLinearSolver_, "cg",
                        Pattern::Selection(linear_solver_options),
                        "linear solvers");

  std::string eigen_solver_options{GetOptionString(kEigenSolverTypeMap_)};
  handler.declare_entry(kEigenSolver_, "pi",
                        Pattern::Selection(eigen_solver_options),
                        "eigenvalue solvers");

  
}

template<typename T>
std::string ParametersDealiiHandler::GetOptionString(
    const std::unordered_map<std::string, T> enum_map) const {
  std::ostringstream oss;
  std::string return_string;
  for (auto const &entry : enum_map) {
    oss << entry.first << "|";
  }
  return_string = oss.str();
  return_string.pop_back();
  return return_string;
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
