#include "parameters_dealii_handler.h"

#include <deal.II/base/parameter_handler.h>

#include <algorithm>
#include <sstream>

namespace bart {

namespace problem {

ParametersDealiiHandler::ParametersDealiiHandler(const std::string filename) {
  dealii::ParameterHandler prm;
  SetUp(prm);
  prm.parse_input(filename, "");
  Parse(prm);
}

void ParametersDealiiHandler::Parse(dealii::ParameterHandler &handler) {

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
  output_filename_base_ = handler.get(key_words_.kOutputFilenameBase_);
  reflective_boundary_ = ParseDealiiMultiple(
      handler.get(key_words_.kReflectiveBoundary_),
      kBoundaryMap_);
  spatial_dimension_ = handler.get_integer(key_words_.kSpatialDimension_);
  spatial_max = ParseDealiiList(handler.get(key_words_.kSpatialMax_));
  transport_model_ = kEquationTypeMap_.at(
      handler.get(key_words_.kTransportModel_));

  // Mesh parameters
  is_mesh_generated_ = handler.get_bool(key_words_.kMeshGenerated_);
  mesh_file_name_ = handler.get(key_words_.kMeshFilename_);
  uniform_refinements_ = handler.get_integer(key_words_.kUniformRefinements_);
  fuel_pin_radius_ = handler.get_double(key_words_.kFuelPinRadius_);
  fuel_pin_triangulation_ = kFuelPinTriangulationTypeMap_.at(
      handler.get(key_words_.kFuelPinTriangulation_));
  is_mesh_pin_resolved_ = handler.get_bool(key_words_.kMeshPinResolved_);

  // Material parameters
  n_materials_ = handler.get_integer(key_words_.kNumberOfMaterials_);
  handler.enter_subsection(key_words_.kMaterialSubsection_);
  material_filenames_ = ParseMap(handler.get(key_words_.kMaterialFilenames_));
  material_map_filename_ = handler.get(key_words_.kMaterialMapFilename_);
  fuel_pin_material_map_filename_ = handler.get(
      key_words_.kFuelPinMaterialMapFilename_);
  handler.leave_subsection();
  
  // Acceleration
  preconditioner_ = kPreconditionerTypeMap_.at(
      handler.get(key_words_.kPreconditioner_));
  block_ssor_factor_ = handler.get_double(key_words_.kBSSOR_Factor_);
  do_nda_ = handler.get_bool(key_words_.kDoNDA_);
  nda_discretization_ = kDiscretizationTypeMap_.at(
      handler.get(key_words_.kNDA_Discretization_));
  nda_linear_solver_ = kLinearSolverTypeMap_.at(
      handler.get(key_words_.kNDALinearSolver_));
  nda_preconditioner_ = kPreconditionerTypeMap_.at(
      handler.get(key_words_.kNDAPreconditioner_));
  nda_block_ssor_factor_ = handler.get_double(key_words_.kNDA_BSSOR_Factor_);
  
  // Solvers
  eigen_solver_ = kEigenSolverTypeMap_.at(
      handler.get(key_words_.kEigenSolver_));
  in_group_solver_ = kInGroupSolverTypeMap_.at(
      handler.get(key_words_.kInGroupSolver_));
  linear_solver_ = kLinearSolverTypeMap_.at(
      handler.get(key_words_.kLinearSolver_));
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
  SetUpMeshParameters(handler);
  SetUpMaterialParameters(handler);
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
  
  handler.declare_entry(key_words_.kOutputFilenameBase_,
                        "bart_output",Pattern::Anything(),
                         "name base of the output file");

  handler.declare_entry(key_words_.kReflectiveBoundary_, "",
                        Pattern::MultipleSelection(
                            GetOptionString(kBoundaryMap_)),
                        "lower case boundary names (xmin, ymax) etc)");
  
  handler.declare_entry(key_words_.kSpatialDimension_, "2",
                        Pattern::Integer(1, 3), "");
  
  try {
    handler.declare_entry(key_words_.kSpatialMax_, "",
                          Pattern::List(Pattern::Double(), 1, 3),
                          "xmax, ymax, zmax of the boundaries, mins are zero");
  } catch (const dealii::ParameterHandler::ExcValueDoesNotMatchPattern &e) {}

  handler.declare_entry(key_words_.kTransportModel_, "none",
                        Pattern::Selection(GetOptionString(kEquationTypeMap_)),
                        "valid names such as ep");
}

// MESH PARAMETERS =============================================================

void ParametersDealiiHandler::SetUpMeshParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;

  handler.declare_entry(key_words_.kMeshGenerated_, "true", Pattern::Bool(),
                        "Boolean to determine if generating mesh in dealii or read in mesh");
  
  handler.declare_entry(key_words_.kMeshFilename_, "", Pattern::Anything(),
                        "mesh file name");
  
  handler.declare_entry(key_words_.kUniformRefinements_, "0",
                        Pattern::Integer(0),
                        "number of uniform refinements desired");

  handler.declare_entry(key_words_.kFuelPinRadius_, "0.5", Pattern::Double(0),
                        "radius of fuel Pin");

  handler.declare_entry(key_words_.kFuelPinTriangulation_, "none",
                        Pattern::Selection(
                            GetOptionString(kFuelPinTriangulationTypeMap_)),
                        "fuel Pin triangulation type");
  handler.declare_entry(key_words_.kMeshPinResolved_, "false", Pattern::Bool(),
                        "Boolean to determine if pPinucing pin-resolved mesh");
                        
}

// MATERIAL PARAMETERS =============================================================

void ParametersDealiiHandler::SetUpMaterialParameters(
    dealii::ParameterHandler &handler) {
  namespace Pattern = dealii::Patterns;

  handler.declare_entry(key_words_.kNumberOfMaterials_, "1", Pattern::Integer(0),
                         "number of materials in the problem");
  
  handler.enter_subsection(key_words_.kMaterialSubsection_);
  
  handler.declare_entry(key_words_.kMaterialMapFilename_, "", Pattern::Anything(),
                        "file name for material id map");
  handler.declare_entry(key_words_.kMaterialFilenames_, "",
                        Pattern::Map(Pattern::Integer(), Pattern::Anything()));
  handler.declare_entry(key_words_.kFuelPinMaterialMapFilename_, "",
                        Pattern::Anything(), "file name for pin material map");
  
  handler.leave_subsection();
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

  handler.declare_entry(key_words_.kNDALinearSolver_, "none",
                        Pattern::Selection(
                            GetOptionString(
                                kLinearSolverTypeMap_,
                                LinearSolverType::kConjugateGradient)),
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
  
  handler.declare_entry(key_words_.kEigenSolver_, "pi",
                        Pattern::Selection(
                            GetOptionString(kEigenSolverTypeMap_)),
                        "eigenvalue solvers");

  handler.declare_entry(key_words_.kInGroupSolver_, "si",
                        Pattern::Selection(
                            GetOptionString(kInGroupSolverTypeMap_)),
                        "in-group solvers");
  
  handler.declare_entry(key_words_.kLinearSolver_, "cg",
                        Pattern::Selection(
                            GetOptionString(kLinearSolverTypeMap_)),
                        "linear solvers");
  
  handler.declare_entry(key_words_.kMultiGroupSolver_, "gs",
                        Pattern::Selection(
                            GetOptionString(kMultiGroupSolverTypeMap_)),
                        "Multi-group solvers");
  
}

void ParametersDealiiHandler::SetUpAngularQuadratureParameters(
    dealii::ParameterHandler &handler) {

  namespace Pattern = dealii::Patterns;
  
  handler.declare_entry(key_words_.kAngularQuad_, "none",
                        Pattern::Selection(
                            GetOptionString(kAngularQuadTypeMap_)),
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

std::unordered_map<int, std::string> ParametersDealiiHandler::ParseMap(
    std::string to_parse) {

  std::unordered_map<int, std::string> return_map;
  
  const std::vector<std::string> pair_strings =
      dealii::Utilities::split_string_list(to_parse, ",");

  for (const std::string& pair_string : pair_strings) {
    const std::vector<std::string> split_pair =
        dealii::Utilities::split_string_list(pair_string, ':');
    return_map[dealii::Utilities::string_to_int(split_pair[0])] = split_pair[1];
  }

  return return_map;
}

} // namespace problem

} // namespace bart
