#include "parameters_dealii_handler.h"

#include <deal.II/base/parameter_handler.h>

#include <algorithm>
#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler() {}

void ParametersDealiiHandler::Parse(const dealii::ParameterHandler &handler) {

  // Parse parameters
  // Basic Parameters
  discretization_ = kDiscretizationTypeMap_.at(
      handler.get(key_words_.kDiscretization_));
  is_eigenvalue_problem_ = handler.get_bool(key_words_.kEigenvalueProblem_);
  have_reflective_bc_ = handler.get_bool(key_words_.kHaveReflectiveBC_);
  fe_polynomial_degree_ = handler.get_integer(key_words_.kFEPolynomialDegree_);
  first_thermal_group_ = handler.get_integer(key_words_.kFirstThermalGroup_);
  n_cells_ = ParseDealiiIntList(handler.get(key_words_.kNCells_));
  n_groups_ = handler.get_integer(key_words_.kNEnergyGroups_);
  n_materials_ = handler.get_integer(key_words_.kNumberOfMaterials_);
  output_filename_base_ = handler.get(key_words_.kOutputFilenameBase_);
  reflective_boundary_ = ParseDealiiMultiple(
      handler.get(key_words_.kReflectiveBoundary_),
      kBoundaryMap_);
  spatial_dimension_ = handler.get_integer(key_words_.kSpatialDimension_);
  spatial_max = ParseDealiiList(handler.get(key_words_.kSpatialMax_));
  transport_model_ = kEquationTypeMap_.at(
      handler.get(key_words_.kTransportModel_));

  // Acceleration
  preconditioner_ = kPreconditionerTypeMap_.at(
      handler.get(key_words_.kPreconditioner_));
  block_ssor_factor_ = handler.get_double(key_words_.kBSSOR_Factor_);
  do_nda_ = handler.get_bool(key_words_.kDoNDA_);
  nda_discretization_ =
      kDiscretizationTypeMap_.at(handler.get(key_words_.kNDA_Discretization_));
  nda_linear_solver_ =
      kLinearSolverTypeMap_.at(handler.get(key_words_.kNDALinearSolver_));
  nda_preconditioner_ = kPreconditionerTypeMap_.at(
      handler.get(key_words_.kNDAPreconditioner_));
  nda_block_ssor_factor_ = handler.get_double(key_words_.kNDA_BSSOR_Factor_);
  
  // Solvers
  eigen_solver_ = kEigenSolverTypeMap_.at(handler.get(key_words_.kEigenSolver_));
  in_group_solver_ = kInGroupSolverTypeMap_.at(handler.get(key_words_.kInGroupSolver_));
  linear_solver_ = kLinearSolverTypeMap_.at(handler.get(key_words_.kLinearSolver_));
  multi_group_solver_ =
      kMultiGroupSolverTypeMap_.at(handler.get(key_words_.kMultiGroupSolver_));

  // Angular Quadrature parameters
  angular_quad_ = kAngularQuadTypeMap_.at(handler.get(key_words_.kAngularQuad_));
  angular_quad_order_ = handler.get_integer(key_words_.kAngularQuadOrder_);
}

void ParametersDealiiHandler::SetUp(dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  /* Set parameters. A try and catch wrapper is required for some. For these,
   * we are using the build-in validity checking for the parameter, but the
   * default values are not valid themselves. Mostly this means that there is
   * no default value.
   */
  SetUpBasicParameters(handler);
  SetUpAccelerationParameters(handler);
  SetUpSolverParameters(handler);
  SetUpAngularQuadratureParameters(handler);
}

/* =============================================================================
 * PARAMETER SET UP
 *
 * This section declares all the entries in the parameter handler that we expect
 * to see in our input file. We use the validation checking that comes with
 * dealii; sometimes default parameters do not meet the validation (if no default
 * parameter makes sense) so a try/catch is needed to ignore the thrown error.
 * The entry is still made. SetUp is divided up to make it a little easier, but
 * it's still a lot of hard to read code.
 * ===========================================================================*/  

void ParametersDealiiHandler::SetUpBasicParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;

  handler.declare_entry(key_words_.kDiscretization_, "cfem",
                        Pattern::Selection(
                            GetOptionString(kDiscretizationTypeMap_)),
                        "HO equation spatial discretization");

  handler.declare_entry(key_words_.kEigenvalueProblem_, "false", Pattern::Bool(),
                        "is problem an eigenvalue problem");

  handler.declare_entry(key_words_.kFEPolynomialDegree_, "1",
                        Pattern::Integer(1),
                        "polynomial degree p for finite element");
  
  handler.declare_entry ("have reflective boundary", "false", Pattern::Bool(),
                         "Does the problem have reflective boundaries");        
  
  handler.declare_entry(key_words_.kFirstThermalGroup_, "0", Pattern::Integer(0),
                        "group number for the first thermal group");
  
  try {
    handler.declare_entry(key_words_.kNCells_ , "",
                          Pattern::List (Pattern::Integer (0), 1, 3),
                          "Geometry is hyper rectangle defined by how many cells exist per direction");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}

  handler.declare_entry(key_words_.kNEnergyGroups_, "1", Pattern::Integer(0),
                        "number of energy groups in the problem");
  
  handler.declare_entry(key_words_.kNumberOfMaterials_, "1", Pattern::Integer(0),
                         "number of materials in the problem");
  
  handler.declare_entry(key_words_.kOutputFilenameBase_, "bart_output",Pattern::Anything(),
                         "name base of the output file");

  std::string boundary_options{GetOptionString(kBoundaryMap_)};
  handler.declare_entry(key_words_.kReflectiveBoundary_, "",
                        Pattern::MultipleSelection(boundary_options),
                        "lower case boundary names (xmin, ymax) etc)");
  
  handler.declare_entry(key_words_.kSpatialDimension_, "2", Pattern::Integer(1, 3), "");
  
  try {
    handler.declare_entry(key_words_.kSpatialMax_, "",
                          Pattern::List(Pattern::Double(), 1, 3),
                          "xmax, ymax, zmax of the boundaries, mins are zero");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}

  std::string equation_options{GetOptionString(kEquationTypeMap_)};
  handler.declare_entry(key_words_.kTransportModel_, "none",
                        Pattern::Selection(equation_options),
                        "valid names such as ep");
}

// ACCELERATION PARAMETERS =====================================================

void ParametersDealiiHandler::SetUpAccelerationParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;

  std::string preconditioner_options{GetOptionString(kPreconditionerTypeMap_)};
  
  handler.declare_entry(key_words_.kPreconditioner_, "amg",
                        Pattern::Selection(preconditioner_options),
                        "Preconditioner");

  handler.declare_entry(key_words_.kBSSOR_Factor_, "1.0", Pattern::Double(0),
                        "damping factor of block SSOR");
  
  handler.declare_entry(key_words_.kDoNDA_, "false", Pattern::Bool(),
                        "Boolean to determine NDA or not");

  handler.declare_entry(key_words_.kNDA_Discretization_, "cfem",
                        Pattern::Selection(
                            GetOptionString(kDiscretizationTypeMap_)),
                        "NDA equation spatial discretization");
  
  // Remove Conjugate Gradient from options for NDA linear solver
  std::string nda_linear_solver_options{GetOptionString(
      kLinearSolverTypeMap_,
      LinearSolverType::kConjugateGradient),};

  handler.declare_entry(key_words_.kNDALinearSolver_, "none",
                        Pattern::Selection(nda_linear_solver_options),
                        "NDA linear solver");

  handler.declare_entry(key_words_.kNDAPreconditioner_, "jacobi",
                        Pattern::Selection(preconditioner_options),
                        "NDA Preconditioner");

  handler.declare_entry(key_words_.kNDA_BSSOR_Factor_, "1.0", Pattern::Double(0),
                        "damping factor of NDA block SSOR");
}

// SOLVER PARAMETERS ===========================================================

void ParametersDealiiHandler::SetUpSolverParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;
  
  std::string eigen_solver_options{GetOptionString(kEigenSolverTypeMap_)};
  handler.declare_entry(key_words_.kEigenSolver_, "pi",
                        Pattern::Selection(eigen_solver_options),
                        "eigenvalue solvers");

  std::string in_group_solver_options{GetOptionString(kInGroupSolverTypeMap_)};
  handler.declare_entry(key_words_.kInGroupSolver_, "si",
                        Pattern::Selection(in_group_solver_options),
                        "in-group solvers");
  
  std::string linear_solver_options{GetOptionString(kLinearSolverTypeMap_)};
  handler.declare_entry(key_words_.kLinearSolver_, "cg",
                        Pattern::Selection(linear_solver_options),
                        "linear solvers");
  
  std::string multigroup_solver_options{GetOptionString(kMultiGroupSolverTypeMap_)};
  handler.declare_entry(key_words_.kMultiGroupSolver_, "gs",
                        Pattern::Selection(multigroup_solver_options),
                        "Multi-group solvers");
  
}

void ParametersDealiiHandler::SetUpAngularQuadratureParameters(
    dealii::ParameterHandler &handler) {

  namespace Pattern = dealii::Patterns;
  
  std::string angular_quad_options{GetOptionString(kAngularQuadTypeMap_)};
  handler.declare_entry(key_words_.kAngularQuad_, "none",
                        Pattern::Selection (angular_quad_options),
                        "angular quadrature types. only LS-GC for multi-D and GL for 1D implemented for now.");
  handler.declare_entry(key_words_.kAngularQuadOrder_, "4", Pattern::Integer(),
                        "Gauss-Chebyshev level-symmetric-like quadrature");
}

template<typename Key>
std::map<Key, bool> ParametersDealiiHandler::ParseDealiiMultiple(
    const std::string to_parse,
    const std::unordered_map<std::string, Key> enum_map) const {

  std::map<Key, bool> return_map;
  
  for (auto const entry : enum_map)
    return_map[entry.second] = false;    

  std::stringstream iss{to_parse};
  std::string read_value;
  
  while(iss >> read_value) {
    // strip commas
    read_value.erase(std::remove(read_value.begin(), read_value.end(), ','),
                     read_value.end());
    
    auto const it = enum_map.find(read_value);
    
    if (it != enum_map.cend())
      return_map[enum_map.at(read_value)] = true;
    
    // Ignore next if comma or space
    if (iss.peek() == ' ')
      iss.ignore();
  }

  return return_map;
  
}

template<typename T>
std::string ParametersDealiiHandler::GetOptionString(
    const std::unordered_map<std::string, T> enum_map,
    const std::vector<T> to_ignore) const {
      
  std::ostringstream oss;
  std::string return_string;
  for (auto const &entry : enum_map) {
    if (std::find(to_ignore.begin(), to_ignore.end(), entry.second) ==
        to_ignore.end())
      oss << entry.first << "|";
  }
  return_string = oss.str();
  if (return_string.size() > 0)
    return_string.pop_back();
  return return_string;
}

template<typename T>
std::string ParametersDealiiHandler::GetOptionString(
    const std::unordered_map<std::string, T> enum_map) const {
  std::vector<T> to_ignore;
  return GetOptionString(enum_map, to_ignore);
}

template<typename T>
std::string ParametersDealiiHandler::GetOptionString(
    const std::unordered_map<std::string, T> enum_map,
    const T to_ignore) const {
  std::vector<T> ignore_vector{to_ignore};
  return GetOptionString(enum_map, ignore_vector);
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
