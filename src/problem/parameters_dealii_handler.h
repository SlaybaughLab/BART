#ifndef BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
#define BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_

#include <string>
#include <unordered_map>

#include <deal.II/base/parameter_handler.h>

#include "parameters_i.h"
#include "parameter_types.h"

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


  // Functions to get problem parameters

  /*! Get linear solver type */
  LinearSolverType    LinearSolver() const override { return linear_solver_; }
  
  /*! Get problem transport model */
  EquationType        TransportModel() const override {
    return transport_model_; }
  
  /*! Get problem spatial dimension */
  std::vector<int>    NCells() const override { return n_cells_; }
  
  /*! Get problem output filename base */
  std::string         OutputFilenameBase() const override {
    return output_filename_base_; }
  
  /*! Get number of spatial dimensions */
  int                 SpatialDimension() const override {
    return spatial_dimension_; }
  
  /*! Get maximum x, y, z size */
  std::vector<double> SpatialMax() const override { return spatial_max; }


 private:
  // Basic parameters  
  EquationType        transport_model_;
  std::vector<int>    n_cells_;
  std::string         output_filename_base_;
  int                 spatial_dimension_;
  std::vector<double> spatial_max;

  // Solvers
  LinearSolverType    linear_solver_;

  // Key-words for input file
  // Basic parameters
  const std::string kNCells_ = "number of cells for x, y, z directions";
  const std::string kOutputFilenameBase_ = "output file name base";
  const std::string kSpatialDimension_ = "problem dimension";
  const std::string kSpatialMax_ = "x, y, z max values of boundary locations";
  const std::string kTransportModel_ = "transport model";

  // Solvers
  const std::string kLinearSolver_ = "ho linear solver name";

  // Options mapping
  const std::unordered_map<std::string, EquationType> kEquationTypeMap_ {
    {"ep",   EquationType::kEvenParity},
    {"saaf", EquationType::kSelfAdjointAngularFlux},
    {"none", EquationType::kNone}
  };

  const std::unordered_map<std::string, LinearSolverType> kLinearSolverTypeMap_ {
    {"cg",    LinearSolverType::kConjugateGradient},
    {"gmres", LinearSolverType::kGMRES},
    {"bicgstab", LinearSolverType::kBiCGSTAB},
    {"direct", LinearSolverType::kDirect}
  };

  // Setup functions
  /*! Set up basic problem parameters */
  void SetUpBasicParameters(dealii::ParameterHandler &handler);
  /*! Set up solver parameters */
  void SetUpSolverParameters(dealii::ParameterHandler &handler);
  
  /*! Parses a ParameterHandler entry of type dealii::Patterns::List with doubles
   * into a vector.
   */
  std::vector<double> ParseDealiiList(std::string to_parse);
  std::vector<int>    ParseDealiiIntList(std::string to_parse);

  /*! Returns a string formed by combining the key strings in a mapping,
   * separated by `|`. Used to generate valid option strings for ParameterHandler
   * entries.
   */
  template<typename T>
  std::string GetOptionString(
      const std::unordered_map<std::string, T> enum_map) const;  

};

} // namespace problem

} // namespace bart

#endif // BART_SRC_PROBLEM_PARAMETERS_DEALII_HANDLER_H_
